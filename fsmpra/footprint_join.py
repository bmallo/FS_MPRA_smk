"""Join existing FiberHMM annotations to a complete molecule ledger by read ID.

All ledger rows survive, including uncalled and unassigned molecules. Local
SQLite bounds memory and enforces one-to-one identities independently of order.
"""
import argparse
import csv
import gzip
import hashlib
import json
import shutil
import time
from collections import Counter, defaultdict
from pathlib import Path

from .io import atomic_output, open_text, read_fasta, read_json, source_info, write_json, write_table
from .processing import check_sources, config_read, extract_interval
from .storage import scratch_database


EXTRA_FIELDS = [
    "footprint_status", "aligned_region_bases", "region_complete",
    "eligible_ref_intervals0", "tf_ref_calls0", "nuc_ref_calls0", "msp_ref_calls0",
    "tf_small_ref_intervals0", "nucleosome_ref_intervals0", "msp_ref_intervals0",
    "tf_small_fraction", "nucleosome_fraction", "msp_fraction",
    "footprint_analysis_eligible", "footprint_analysis_eligible_with_barcode",
]
TRACKS = ("tf", "nuc", "msp")


def packed(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False)


def merge_intervals(intervals):
    """Union overlapping/adjacent half-open intervals without filling gaps."""
    result = []
    for a, b in sorted(intervals):
        if result and a <= result[-1][1]:
            result[-1][1] = max(b, result[-1][1])
        else:
            result.append([a, b])
    return result


def region_blocks(read, start, end):
    """Return (query_start, query_end, reference_start) clipped to the region.

    Stored BAM SEQ and MA share orientation, including reverse alignments.
    Insertions and clips have no reference bases; deletions remain gaps.
    """
    q, r = 0, read.reference_start
    blocks = []
    for op, size in read.cigartuples or []:
        if op in (0, 7, 8):
            a, b = max(r, start), min(r + size, end)
            if a < b:
                blocks.append((q + a - r, q + b - r, a))
            q += size
            r += size
        elif op in (1, 4):
            q += size
        elif op in (2, 3):
            r += size
        elif op not in (5, 6):
            raise ValueError("Unsupported CIGAR operation")
        if r >= end:
            break
    return blocks


def project_call(start, length, blocks):
    spans = []
    for qa, qb, ra in blocks:
        a, b = max(start, qa), min(start + length, qb)
        if a < b:
            spans.append([ra + a - qa, ra + b - qa])
    return merge_intervals(spans)


def region_calls(read, blocks):
    """Validate MA/AQ/AN and retain every call intersecting the region.

    Supported FiberHMM quality specifications use byte-valued Q fields. Quality
    bytes and optional annotation names are preserved, not used for new filters.
    Absent AQ is explicitly null; absent sparse tracks are evaluated empty.
    """
    pieces = read.get_tag("MA").split(";")
    if int(pieces[0]) != read.query_length:
        raise ValueError("MA/query-length mismatch")
    aq = read.get_tag("AQ") if read.has_tag("AQ") else None
    if aq is not None and (getattr(aq, "typecode", None) != "B"):
        raise ValueError("Expected unsigned byte AQ array")
    names = read.get_tag("AN").split(",") if read.has_tag("AN") else None
    result = {name: [] for name in TRACKS}
    seen, cursor, annotation = set(), 0, 0
    for section in pieces[1:]:
        if not section:
            continue
        header, sep, body = section.partition(":")
        track, plus, quality_spec = header.partition("+")
        if not sep or not plus or track not in TRACKS or set(quality_spec) - {"Q"} or track in seen:
            raise ValueError("Unsupported or repeated FiberHMM MA section: " + header)
        seen.add(track)
        width = len(quality_spec)
        for token in body.split(",") if body else []:
            start, length = map(int, token.split("-"))
            start -= 1
            if start < 0 or length <= 0 or start + length > read.query_length:
                raise ValueError("MA interval outside stored SEQ")
            if aq is not None and cursor + width > len(aq):
                raise ValueError("AQ shorter than MA quality specification")
            if names is not None and annotation >= len(names):
                raise ValueError("AN shorter than MA annotation count")
            spans = project_call(start, length, blocks)
            if spans:
                result[track].append({"query_start0": start, "query_length": length,
                    "ref_intervals0": spans, "quality_spec": quality_spec,
                    "quality": list(aq[cursor:cursor + width]) if aq is not None else ([] if not width else None),
                    "name": names[annotation] if names is not None else None})
            cursor += width
            annotation += 1
    if aq is not None and cursor != len(aq):
        raise ValueError("AQ longer than MA quality specification")
    if names is not None and annotation != len(names):
        raise ValueError("AN longer than MA annotation count")
    return result, seen


def annotate(read, start, end, tf_min=10, tf_max=60, nuc_min=85):
    blocks = region_blocks(read, start, end)
    aligned = merge_intervals([[r, r + qb - qa] for qa, qb, r in blocks])
    n = sum(b - a for a, b in aligned)
    row = dict.fromkeys(EXTRA_FIELDS, "")
    row.update(footprint_status="missing_MA", aligned_region_bases=n,
               region_complete=int(n == end - start), eligible_ref_intervals0=packed(aligned),
               footprint_analysis_eligible=0, footprint_analysis_eligible_with_barcode=0)
    if not read.has_tag("MA"):
        return row, blocks, set()
    calls, present = region_calls(read, blocks)
    row["footprint_status"] = "evaluated"
    for track in TRACKS:
        row[track + "_ref_calls0"] = packed(calls[track])
    for label, track, lo, hi in (("tf_small", "tf", tf_min, tf_max),
                                ("nucleosome", "nuc", nuc_min, None), ("msp", "msp", 1, None)):
        spans = merge_intervals([span for call in calls[track]
            if call["query_length"] >= lo and (hi is None or call["query_length"] <= hi)
            for span in call["ref_intervals0"]])
        row[label + "_ref_intervals0"] = packed(spans)
        if n == end - start:
            row[label + "_fraction"] = sum(b - a for a, b in spans) / n
    return row, blocks, present


def join_footprints(config, catalog, ledger, bam, contract, output,
                    cache_mib=128, max_disk_gib=40, tf_min=10, tf_max=60, nuc_min=85):
    import pysam
    if not 1 <= tf_min <= tf_max or nuc_min < 1:
        raise ValueError("Invalid footprint size thresholds")
    cfg, lib, spec = config_read(config), read_json(catalog)["library"], read_json(contract)
    if spec.get("format") != "fiberhmm_ma_sparse" or set(spec.get("evaluated_tracks", [])) != set(TRACKS):
        raise ValueError("An explicit FiberHMM sparse evaluated-track contract is required")
    if Path(spec["bam"]).resolve() != Path(bam).resolve() or spec["sample_id"] != cfg["sample"]["id"]:
        raise ValueError("Footprint contract does not identify this sample/BAM")
    ledger = Path(ledger)
    validation_path = ledger.parent / "validation.json"
    validation = read_json(validation_path)
    if validation.get("schema_version") != "fsmpra.ledger.v1":
        raise ValueError("Expected a validated complete ledger")
    if cfg["library"]["id"] != lib["id"]:
        raise ValueError("Configuration/catalog library mismatch")
    for path in (config, catalog):
        actual = source_info(path, True)
        recorded = next((s for s in validation["sources"] if s["path"] == actual["path"]), None)
        if recorded != actual:
            raise ValueError("Ledger provenance disagrees with config/catalog")
    start, end = lib["mutable_ref_start0"], lib["mutable_ref_end0"]
    reference_path = Path(cfg["library"]["reference_fasta"])
    if not reference_path.is_absolute():
        reference_path = Path(cfg["project_root"]) / reference_path
    reference = read_fasta(reference_path, lib["contig"])
    if hashlib.sha256(reference.encode()).hexdigest() != lib["reference_sequence_sha256"]:
        raise ValueError("Reference sequence disagrees with catalog")
    variants = {v["variant_id"]: v for v in lib["variants"]}
    sources = [source_info(p, True) for p in (config, catalog, contract, validation_path, reference_path,
               __file__, Path(__file__).with_name("processing.py"), Path(__file__).with_name("storage.py"))]
    sources += [source_info(ledger), source_info(bam), source_info(str(bam) + ".bai")]
    evidence = spec.get("evidence", {})
    sources += [source_info(p, True) for p in evidence.get("caller_logs", [])]
    if evidence.get("sparse_writer"):
        sources.append(source_info(evidence["sparse_writer"], True))
    if shutil.disk_usage(cfg["scratch_dir"]).free < (max_disk_gib + 2) * 1024**3:
        raise ValueError("Insufficient local scratch for requested database cap")
    counts, by_variant = Counter(), defaultdict(Counter)
    started = time.monotonic()
    with atomic_output(output) as out, scratch_database(cfg["scratch_dir"], cache_mib, max_disk_gib) as db:
        db.execute("CREATE TABLE ledger (id TEXT PRIMARY KEY, row TEXT NOT NULL, seen INTEGER NOT NULL DEFAULT 0)")
        with open_text(ledger) as src:
            reader = csv.reader(src, delimiter="\t")
            fields = next(reader)
            required = {"read_id", "sample_id", "library_id", "observed_sequence", "observed_variant_id", "barcode",
                        "read_length", "strand", "analysis_eligible", "assigned_variant_id", "corrected_barcode"}
            if not required <= set(fields) or len(fields) != len(set(fields)) or set(fields) & set(EXTRA_FIELDS):
                raise ValueError("Invalid ledger columns")
            idx = {name: i for i, name in enumerate(fields)}
            for values in reader:
                if len(values) != len(fields) or any("\t" in x or "\n" in x or "\r" in x for x in values):
                    raise ValueError("Invalid ledger row")
                if values[idx["sample_id"]] != cfg["sample"]["id"] or values[idx["library_id"]] != lib["id"]:
                    raise ValueError("Ledger sample/library mismatch")
                vid, assigned = values[idx["observed_variant_id"]], values[idx["assigned_variant_id"]]
                if (vid and vid not in variants) or (assigned and (assigned != vid or not values[idx["corrected_barcode"]])):
                    raise ValueError("Invalid ledger genotype assignment")
                db.execute("INSERT INTO ledger(id,row) VALUES (?,?)", (values[idx["read_id"]], "\t".join(values)))
                counts["ledger_molecules"] += 1
                if counts["ledger_molecules"] % 10000 == 0:
                    db.commit()
                if counts["ledger_molecules"] % 500000 == 0:
                    print(f"Indexed ledger: {counts['ledger_molecules']:,}", flush=True)
        db.commit()
        if counts["ledger_molecules"] != validation["counts"]["molecules"]:
            raise ValueError("Ledger count differs from validated count")

        def emit(values, annotation, writer):
            writer.writerow(values + [annotation.get(f, "") for f in EXTRA_FIELDS])
            counts["output_molecules"] += 1
            groups = [counts, by_variant[values[idx["observed_variant_id"]]]]
            for group in groups:
                group[annotation["footprint_status"]] += 1
                group["region_complete"] += int(annotation.get("region_complete", 0) or 0)
                group["ledger_analysis_eligible"] += int(values[idx["analysis_eligible"]])
                group["footprint_analysis_eligible"] += annotation["footprint_analysis_eligible"]
                group["footprint_analysis_eligible_with_barcode"] += annotation["footprint_analysis_eligible_with_barcode"]
            by_variant[values[idx["observed_variant_id"]]]["ledger_molecules"] += 1

        with gzip.open(out / "molecule_footprints.tsv.gz", "wt", compresslevel=1, newline="") as dest, \
                pysam.AlignmentFile(bam, "rb", threads=1) as src:
            if src.get_reference_length(lib["contig"]) != lib["reference_length"]:
                raise ValueError("Footprint BAM reference length mismatch")
            indexed = next(x.mapped for x in src.get_index_statistics() if x.contig == lib["contig"])
            programs = src.header.to_dict().get("PG", [])
            writer = csv.writer(dest, delimiter="\t")
            writer.writerow(fields + EXTRA_FIELDS)
            for read in src.fetch(lib["contig"]):
                counts["bam_alignment_records"] += 1
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    counts["bam_nonprimary_excluded"] += 1
                    continue
                found = db.execute("SELECT row,seen FROM ledger WHERE id=?", (read.query_name,)).fetchone()
                if found is None or found[1]:
                    raise ValueError("Unknown or duplicate primary footprint read: " + read.query_name)
                values = found[0].split("\t")
                if int(values[idx["read_length"]]) != read.query_length or values[idx["strand"]] != ("-" if read.is_reverse else "+"):
                    raise ValueError("Ledger/BAM length or orientation mismatch: " + read.query_name)
                sequence = extract_interval(read, start, end, reference)[0]
                barcode = extract_interval(read, lib["barcode_start0"], lib["barcode_end0"], reference)[0]
                if lib["barcode_orientation"] == "reverse_complement":
                    barcode = barcode.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
                if sequence != values[idx["observed_sequence"]] or barcode != values[idx["barcode"]]:
                    raise ValueError("Ledger/BAM mutable sequence or barcode mismatch: " + read.query_name)
                annotation, blocks, present = annotate(read, start, end, tf_min, tf_max, nuc_min)
                if annotation["footprint_status"] == "evaluated":
                    counts["evaluated_without_AQ"] += int(not read.has_tag("AQ"))
                    for track in set(TRACKS) - present:
                        counts["omitted_empty_track_" + track] += 1
                if values[idx["analysis_eligible"]] == "1" and annotation["footprint_status"] == "evaluated" and annotation["region_complete"]:
                    variant = variants[values[idx["observed_variant_id"]]]
                    projected_seq = "".join(read.query_sequence[a:b] for a, b, _ in blocks)
                    if variant["class"] not in ("WT", "SNV") or sequence != variant["sequence"]:
                        raise ValueError("Analysis-eligible ledger genotype disagrees with catalog")
                    if projected_seq == variant["sequence"]:
                        annotation["footprint_analysis_eligible"] = 1
                        annotation["footprint_analysis_eligible_with_barcode"] = int(bool(values[idx["assigned_variant_id"]]))
                    else:
                        counts["eligible_projection_sequence_disagreement"] += 1
                emit(values, annotation, writer)
                db.execute("UPDATE ledger SET seen=1 WHERE id=?", (read.query_name,))
                counts["matched_primary_molecules"] += 1
                if counts["matched_primary_molecules"] % 10000 == 0:
                    db.commit()
                if counts["matched_primary_molecules"] % 250000 == 0:
                    print(f"Footprints joined: {counts['matched_primary_molecules']:,}", flush=True)
            if counts["bam_alignment_records"] != indexed:
                raise ValueError("Footprint traversal differs from BAM index")
            db.commit()
            missing = dict.fromkeys(EXTRA_FIELDS, "")
            missing.update(footprint_status="missing_bam", footprint_analysis_eligible=0,
                           footprint_analysis_eligible_with_barcode=0)
            for (row,) in db.execute("SELECT row FROM ledger WHERE seen=0"):
                emit(row.split("\t"), missing, writer)
        if counts["output_molecules"] != counts["ledger_molecules"]:
            raise ValueError("Footprint join lost ledger molecules")
        if counts["ledger_analysis_eligible"] != validation["counts"]["analysis_eligible"]:
            raise ValueError("Analysis cohort count changed")
        coverage_fields = ["variant_id", "source_variant", "genotype_class", "ledger_molecules", "evaluated", "missing_MA",
            "missing_bam", "region_complete", "ledger_analysis_eligible", "footprint_analysis_eligible",
            "footprint_analysis_eligible_with_barcode"]
        rows = []
        for vid in list(variants) + [""]:
            variant = variants.get(vid, {"source_variant": "unresolved", "class": "unresolved"})
            rows.append(dict(variant_id=vid, source_variant=variant["source_variant"], genotype_class=variant["class"],
                             **{k: by_variant[vid][k] for k in coverage_fields[3:]}))
        write_table(out / "variant_footprint_coverage.tsv", rows, coverage_fields)
        check_sources(sources)
        report = {"schema_version": "fsmpra.footprint_join.v1", "sources": sources, "counts": dict(counts),
            "sample_id": cfg["sample"]["id"], "library_id": lib["id"], "interval": [start, end],
            "coordinates": "Reference coordinates are 0-based half-open; query_start0 uses stored BAM SEQ orientation.",
            "footprint_contract": spec, "programs": programs, "seconds": time.monotonic() - started,
            "size_thresholds": {"tf_min": tf_min, "tf_max": tf_max, "nuc_min": nuc_min},
            "storage": {"cache_mib": cache_mib, "max_disk_gib": max_disk_gib, "scratch_dir": cfg["scratch_dir"]},
            "missingness": "Absent BAM/MA is missing, never zero. Omitted tracks in present MA are evaluated empty under the explicit contract.",
            "call_fields": "Each region-overlapping call retains query_start0, original query_length, clipped ref_intervals0, quality_spec, raw AQ quality bytes (null if absent), and AN name (null if absent).",
            "limitations": ["Existing caller output reused; no new footprint calling or occupancy inference.",
                "Fractions require all reference bases aligned; eligibility additionally requires exact WT/SNV projected sequence.",
                "Size filters use MA piece length; circular split pieces retain names when present but are not reunited.",
                "Barcode assignments are preserved; molecules are not collapsed by barcode."]}
        write_json(out / "audit.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("config", "catalog", "ledger", "bam", "contract", "output"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--cache-mib", type=int, default=128)
    parser.add_argument("--max-disk-gib", type=float, default=40)
    parser.add_argument("--tf-min", type=int, default=10)
    parser.add_argument("--tf-max", type=int, default=60)
    parser.add_argument("--nuc-min", type=int, default=85)
    args = vars(parser.parse_args())
    print(json.dumps(join_footprints(**args)["counts"], indent=2))
