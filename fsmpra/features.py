"""Descriptive, allele-linked FiberHMM measurements for the initial pilot.

MSP coverage is a caller-derived accessibility proxy, not an absolute probability
of biological accessibility. No hypothesis tests or calibrated calls are made.
"""

import csv
import hashlib
import logging
from collections import Counter
from pathlib import Path

from . import SCHEMA_VERSION, __version__
from .io import new_output, open_text, read_json, source_info, write_json, write_table


def parse_ma(value, query_length):
    """Strict FiberHMM MA parser; starts are 1-based stored-SEQ coordinates."""
    parts = value.split(";")
    if int(parts[0]) != query_length:
        raise ValueError("MA/query-length mismatch")
    tracks = {}
    for segment in parts[1:]:
        if not segment:
            continue
        label, sep, body = segment.partition(":")
        if not sep:
            raise ValueError("Malformed MA section")
        track = label.split("+", 1)[0]
        if track in tracks:
            raise ValueError(f"Duplicate MA track {track}")
        intervals = []
        for item in body.split(",") if body else []:
            start, length = map(int, item.split("-"))
            start -= 1
            if start < 0 or length <= 0 or start + length > query_length:
                raise ValueError("MA interval outside query")
            intervals.append((start, length))
        tracks[track] = intervals
    return tracks


def projected_coverage(intervals, qmap, start0, end0, min_length=1, max_length=None):
    """Cover only aligned query bases; never fill reference-deletion gaps."""
    covered = set()
    lengths = []
    for start, length in intervals:
        if length < min_length or (max_length is not None and length > max_length):
            continue
        positions = {r - start0 for r in qmap[start:start + length]
                     if r is not None and start0 <= r < end0}
        if positions:
            lengths.append(length)
            covered.update(positions)
    return covered, lengths


def interval_runs(positions):
    """Zero-based half-open runs; input positions are relative to the region."""
    result = []
    for p in sorted(positions):
        if result and p == result[-1][1]:
            result[-1][1] = p + 1
        else:
            result.append([p, p + 1])
    return result


def divide(numerators, denominators):
    return [round(a / b, 6) if b else None for a, b in zip(numerators, denominators)]


def masked_mean(counts, n_reads, positions):
    return sum(counts[i] for i in positions) / (n_reads * len(positions)) if n_reads and positions else None


def measure(catalog_path, molecules_dir, output, min_passes=10, tf_min=10, tf_max=60,
            nuc_min=85, mask_radius=2, examples_per_variant=4):
    import pysam

    if min_passes < 0 or not 1 <= tf_min <= tf_max or nuc_min < 1 or mask_radius < 0 or examples_per_variant < 0:
        raise ValueError("Invalid measurement parameters")
    if Path(output).exists():
        raise FileExistsError(output)
    catalog = read_json(catalog_path)
    mol_dir = Path(molecules_dir)
    audit = read_json(mol_dir / "audit.json")
    if source_info(catalog_path, True) != audit["catalog"]:
        raise ValueError("Molecule import and catalog provenance disagree")
    sample = next(s for s in catalog["samples"] if s["id"] == audit["sample_id"])
    contract = sample.get("footprint_contract", {})
    evaluated_tracks = set(contract.get("evaluated_tracks", []))
    if evaluated_tracks and contract.get("format") != "fiberhmm_ma_sparse":
        raise ValueError("Evaluated-track declaration requires a supported footprint format")
    if evaluated_tracks - {"tf", "msp", "nuc"}:
        raise ValueError("Unknown evaluated track")
    lib = next(l for l in catalog["libraries"] if l["id"] == sample["library_id"])
    start, end = lib["mutable_ref_start0"], lib["mutable_ref_end0"]
    length = end - start
    bam_path = sample["inputs"]["footprints_bam"]["path"]
    bam_source = source_info(bam_path)
    if bam_source != sample["inputs"]["footprints_bam"]:
        raise ValueError("BAM changed since catalog creation")
    molecule_source = source_info(mol_dir / "molecules.tsv.gz", True)
    by_id = {v["variant_id"]: v for v in lib["variants"]}
    specs = {"tf_small": ("tf", tf_min, tf_max), "nucleosome": ("nuc", nuc_min, None),
             "msp": ("msp", 1, None)}
    counts = Counter()
    eligible = {}
    with open_text(mol_dir / "molecules.tsv.gz") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts["imported"] += 1
            # One primary exclusion reason per molecule; no methylation selection.
            reason = ("missing_metrics" if row["metrics_present"] != "1" else
                      "not_full_plasmid" if row["plasmid_class"] != "plasmid_full" else
                      "not_WT_or_single_SNV" if row["genotype_class"] not in ("WT", "SNV") else
                      "low_passes" if int(row["passes"]) < min_passes else None)
            if reason:
                counts[reason] += 1
                continue
            if row["variant_id"] not in by_id:
                raise ValueError("Unknown molecule variant_id")
            if row["read_id"] in eligible:
                raise ValueError("Duplicate molecule ID")
            eligible[row["read_id"]] = row
    logging.info("Measurement cohort: %s sampled WT/SNV molecules", f"{len(eligible):,}")
    aggregate = {}
    for v in lib["variants"]:
        if v["class"] not in ("WT", "SNV"):
            continue
        aggregate[v["variant_id"]] = {
            "n_reads": 0, "n_forward": 0, "n_reverse": 0,
            "methylation_histogram": [0] * 20, "msp_fraction_histogram": [0] * 20,
            "msp_fraction_n": 0, "msp_fraction_sum": 0.0,
            "examples": [], "tracks": {key: {"called_reads": [0] * length, "eligible_reads": [0] * length,
                                             "full_coverage_counts": [0] * length, "n_full_coverage": 0,
                                             "footprint_length_histogram": Counter()}
                                        for key in specs}}
    seen = set()
    feature_rows = []
    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam:
        if bam.get_reference_length(lib["contig"]) != lib["reference_length"]:
            raise ValueError("BAM reference length disagrees with catalog")
        for i, read in enumerate(bam.fetch(lib["contig"], start, end), 1):
            if i % 1000000 == 0:
                logging.info("BAM: %s alignments scanned, %s sampled reads measured", f"{i:,}", f"{counts['measured']:,}")
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            mol = eligible.get(read.query_name)
            if mol is None:
                continue
            if read.query_name in seen:
                raise ValueError(f"Duplicate primary alignment {read.query_name}")
            seen.add(read.query_name)
            qmap = read.get_reference_positions(full_length=True)
            pairs = [(q, r - start) for q, r in enumerate(qmap) if r is not None and start <= r < end]
            aligned = {p for q, p in pairs}
            variant = by_id[mol["variant_id"]]
            # Full-region summaries require an exact designed WT/SNV sequence.
            # This also verifies the ledger/BAM join and protects the coordinate contract.
            if len(aligned) != length:
                counts["incomplete_analysis_interval"] += 1
                continue
            qlo, qhi = pairs[0][0], pairs[-1][0] + 1
            if read.query_sequence[qlo:qhi].upper() != variant["sequence"]:
                counts["sequence_disagrees_with_genotype"] += 1
                continue
            if not read.has_tag("MA"):
                counts["missing_MA"] += 1
                continue
            tracks = parse_ma(read.get_tag("MA"), read.query_length)
            # FiberHMM's sparse writer omits empty sections. Only interpret that
            # as zero for tracks explicitly documented as evaluated for this input.
            for track in evaluated_tracks:
                if track not in tracks:
                    counts[f"omitted_empty_track_{track}"] += 1
                    tracks[track] = []
            item = aggregate[mol["variant_id"]]
            item["n_reads"] += 1
            item["n_reverse" if read.is_reverse else "n_forward"] += 1
            mf = float(mol["global_methylation_fraction"])
            item["methylation_histogram"][min(19, int(mf * 20))] += 1
            count_row = {"sample_id": sample["id"], "read_id": read.query_name,
                         "variant_id": mol["variant_id"], "source_variant": mol["source_variant"],
                         "global_methylation_fraction": mf, "strand": "-" if read.is_reverse else "+",
                         "aligned_region_bases": len(aligned)}
            display_tracks = {}
            for name, (track, lo, hi) in specs.items():
                if track not in tracks:
                    counts[f"missing_track_{track}"] += 1
                    count_row[f"{name}_fraction"] = None
                    count_row[f"{name}_segments"] = None
                    display_tracks[name] = None
                    continue
                covered, lengths = projected_coverage(tracks[track], qmap, start, end, lo, hi)
                group = item["tracks"][name]
                for p in aligned:
                    group["eligible_reads"][p] += 1
                for p in covered:
                    group["called_reads"][p] += 1
                    group["full_coverage_counts"][p] += 1
                group["n_full_coverage"] += 1
                group["footprint_length_histogram"].update(lengths)
                fraction = len(covered) / length
                count_row[f"{name}_fraction"] = fraction
                count_row[f"{name}_segments"] = len(lengths)
                display_tracks[name] = interval_runs(covered)
                if name == "msp":
                    item["msp_fraction_histogram"][min(19, int(fraction * 20))] += 1
                    item["msp_fraction_n"] += 1
                    item["msp_fraction_sum"] += fraction
            feature_rows.append(count_row)
            # Representative examples: deterministic bottom-k by read hash, not BAM order.
            if examples_per_variant:
                rank = hashlib.sha256(read.query_name.encode()).hexdigest()
                examples = item["examples"]
                if len(examples) < examples_per_variant or rank < examples[-1]["rank"]:
                    examples.append({"rank": rank, "read_id": read.query_name, "strand": count_row["strand"],
                                     "global_methylation_fraction": mf, "tracks": display_tracks})
                    examples.sort(key=lambda x: x["rank"])
                    del examples[examples_per_variant:]
            counts["measured"] += 1
    counts["eligible_missing_from_BAM_region"] = len(eligible) - len(seen)
    if source_info(bam_path) != bam_source:
        raise ValueError("BAM changed during measurement")
    wt_id = next(v["variant_id"] for v in lib["variants"] if v["class"] == "WT")
    wt = aggregate[wt_id]
    summary, records = [], []
    for vid, item in aggregate.items():
        variant = by_id[vid]
        offset = variant["ref_pos0"] - start if variant["class"] == "SNV" else None
        common = [p for p in range(length) if offset is None or abs(p - offset) > mask_radius]
        record = {"variant_id": vid, "source_variant": variant["source_variant"], "class": variant["class"],
                  "oligo_pos1": variant["oligo_pos1"], "genomic_pos1": variant["genomic_pos1"],
                  "ref": variant["ref"], "alt": variant["alt"], **item}
        row = {"sample_id": sample["id"], "variant_id": vid, "source_variant": variant["source_variant"],
               "class": variant["class"], "oligo_pos1": variant["oligo_pos1"], "n_reads": item["n_reads"],
               "n_forward": item["n_forward"], "n_reverse": item["n_reverse"],
               "msp_fraction_mean": item["msp_fraction_sum"] / item["msp_fraction_n"] if item["msp_fraction_n"] else None,
               "msp_fraction_n": item["msp_fraction_n"]}
        for name in specs:
            group, baseline = item["tracks"][name], wt["tracks"][name]
            occ = divide(group["called_reads"], group["eligible_reads"])
            wt_occ = divide(baseline["called_reads"], baseline["eligible_reads"])
            group["occupancy"] = occ
            group["delta_vs_WT"] = [round(a - b, 6) if a is not None and b is not None else None for a, b in zip(occ, wt_occ)]
            group["footprint_length_histogram"] = dict(sorted(group["footprint_length_histogram"].items()))
            a = masked_mean(group["full_coverage_counts"], group["n_full_coverage"], common)
            b = masked_mean(baseline["full_coverage_counts"], baseline["n_full_coverage"], common)
            row[f"{name}_common_site_mean"] = a
            row[f"{name}_WT_common_site_mean"] = b
            row[f"{name}_common_site_delta"] = a - b if a is not None and b is not None else None
        record["summary"] = row
        summary.append(row)
        records.append(record)
    result = {"schema_version": SCHEMA_VERSION, "software_version": __version__,
              "sample": {k: v for k, v in sample.items() if k != "inputs"},
              "library_id": lib["id"], "reference_sequence_sha256": lib["reference_sequence_sha256"],
              "region": {"contig": lib["contig"], "start0": start, "end0": end,
                         "oligo_start1": lib["oligo_mutable_start1"], "sequence": lib["mutable_sequence"],
                         "status": lib["analysis_interval_status"]},
              "sampling": audit["sampling"], "counts": dict(counts),
              "method": {"caller": "existing FiberHMM MA", "MA_coordinates": "1-based stored-SEQ; no reverse-strand flip",
                         "footprint_contract": contract,
                         "tracks": specs, "min_passes": min_passes, "methylation_selection": None,
                         "common_site_exclusion_radius_bp": mask_radius, "histogram_edges": [i / 20 for i in range(21)],
                         "inference": "descriptive_only", "cohort": "full plasmid, WT or single SNV, pass threshold, exact designed mutable sequence"},
              "provenance": {"catalog": source_info(catalog_path, True), "molecules": molecule_source, "bam": bam_source},
              "limitations": ["Pilot sample; no p-values, confidence intervals, or claims of discovery.",
                              "MSP coverage is a caller-derived accessibility proxy, not absolute accessibility.",
                              "Sequence-dependent footprint observability is not calibrated; masking is a sensitivity readout, not a validated correction.",
                              "Only the 211-bp mutable region is summarized; promoter boundary selection remains provisional.",
                              "Exact sequence validation is confined to the mutable interval; flanking mutations remain unassessed.",
                              "Absent MA is missing data. Omitted tracks are zero only when the manifest declares they were evaluated by the sparse writer. All methylation levels remain represented."],
              "variants": records,
              "mpra": [r for r in catalog["mpra"] if r["variant_id"] in aggregate]}
    out = new_output(output)
    write_json(out / "pilot.json", result)
    write_table(out / "variant_summary.tsv", summary, list(summary[0]))
    fields = ["sample_id", "read_id", "variant_id", "source_variant", "global_methylation_fraction", "strand", "aligned_region_bases"]
    fields += [field for name in specs for field in (f"{name}_fraction", f"{name}_segments")]
    write_table(out / "molecule_features.tsv.gz", feature_rows, fields)
    write_json(out / "audit.json", {k: result[k] for k in ("schema_version", "sample", "sampling", "counts", "method", "provenance", "limitations")})
    return {"sample_id": sample["id"], "counts": dict(counts), "output": str(out)}
