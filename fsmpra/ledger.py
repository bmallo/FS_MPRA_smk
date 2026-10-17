"""Complete plasmid ledger from an existing indexed alignment.

Bounded local BAM chunks allow extraction on allocated CPUs. Outputs remain
molecule-level; barcode grouping never deduplicates physical reads.
"""
import argparse
import csv
import gzip
import itertools
import multiprocessing
import shutil
import tempfile
import time
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .io import atomic_output, open_text, read_json, source_info, write_json
from .processing import OBS_FIELDS, SCHEMA, check_sources, config_read, extract


def _chunk_extract(config, catalog, bam, output, first):
    return extract(config, catalog, bam, output, write_header=first)


def extract_complete(config, catalog, bam, output, workers=3, allocated_cpus=4, chunk_reads=100000):
    import pysam
    if workers < 1 or workers + 1 > allocated_cpus or chunk_reads < 1:
        raise ValueError("Allocate one CPU for the reader plus one per extraction worker; chunk_reads must be positive")
    cfg = config_read(config)
    lib = read_json(catalog)["library"]
    if shutil.disk_usage(cfg["scratch_dir"]).free < 10 * 1024**3:
        raise ValueError("At least 10 GiB free local scratch required")
    sources = [source_info(config, True), source_info(catalog, True), source_info(bam), source_info(str(bam) + ".bai"),
               source_info(__file__, True), source_info(Path(__file__).with_name("processing.py"), True)]
    counts, dispatch = Counter(), Counter()
    started = time.monotonic()
    chunk_audit = []
    with atomic_output(output) as out, tempfile.TemporaryDirectory(prefix="fsmpra-ledger-", dir=cfg["scratch_dir"]) as temp:
        with pysam.AlignmentFile(bam, "rb", threads=1) as src:
            if src.get_reference_length(lib["contig"]) != lib["reference_length"]:
                raise ValueError("Alignment reference length disagrees with library")
            indexed = next(x for x in src.get_index_statistics() if x.contig == lib["contig"])
            programs = src.header.to_dict().get("PG", [])
            pending = deque()
            # Concatenated gzip members are valid gzip. Only the first shard
            # has a TSV header, so downstream readers see one ordinary table.
            with open(out / "observations.tsv.gz", "wb") as combined, ProcessPoolExecutor(
                    max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:

                def collect():
                    future, path, chunk_dir, index, n = pending.popleft()
                    report = future.result()
                    if report["counts"].get("plasmid_primary", 0) != n:
                        raise ValueError("Chunk molecule count changed")
                    with open(chunk_dir / "observations.tsv.gz", "rb") as handle:
                        shutil.copyfileobj(handle, combined, length=1024*1024)
                    counts.update(report["counts"])
                    chunk_audit.append({"chunk": index, "molecules": n, "counts": report["counts"]})
                    path.unlink()
                    shutil.rmtree(chunk_dir)
                    print(f"Extracted {counts['plasmid_primary']:,} primary plasmids ({len(chunk_audit)} chunks)", flush=True)

                target = src.fetch(lib["contig"])
                index, done = 0, False
                while not done:
                    path = Path(temp) / f"chunk-{index:05d}.bam"
                    n = 0
                    with pysam.AlignmentFile(path, "wb0", template=src, threads=1) as dest:
                        while n < chunk_reads:
                            read = next(target, None)
                            if read is None:
                                done = True
                                break
                            dispatch["target_alignment_records"] += 1
                            if read.is_secondary or read.is_supplementary:
                                dispatch["secondary_or_supplementary"] += 1
                                continue
                            if read.is_unmapped:
                                raise ValueError("Unexpected unmapped record in indexed contig traversal")
                            dest.write(read)
                            n += 1
                    if n:
                        chunk_dir = Path(temp) / f"extracted-{index:05d}"
                        future = pool.submit(_chunk_extract, str(config), str(catalog), str(path), str(chunk_dir), index == 0)
                        pending.append((future, path, chunk_dir, index, n))
                        index += 1
                        if len(pending) >= workers:
                            collect()
                    else:
                        path.unlink()
                while pending:
                    collect()
            if not chunk_audit:
                raise ValueError("No primary plasmid alignments")
            if dispatch["target_alignment_records"] != indexed.mapped:
                raise ValueError("Traversal count disagrees with BAM index")
        counts["alignment_records"] = dispatch["target_alignment_records"]
        counts["secondary_or_supplementary"] = dispatch["secondary_or_supplementary"]
        check_sources(sources)
        report = {"schema_version": SCHEMA, "sources": sources, "counts": dict(counts), "qc": cfg["qc"],
            "table_header": True, "chunks": chunk_audit, "workers": workers, "allocated_cpus": allocated_cpus,
            "seconds": time.monotonic() - started, "programs": programs,
            "scope": "All primary alignments on the configured plasmid contig; other genomic contigs/unmapped reads were not traversed.",
            "barcode_orientation": lib["barcode_orientation"],
            "genotype_method": "Exact observed mutable sequence; no barcode-based genotype rescue",
            "global_duplicate_check": "Barcode ingestion enforces a primary key across the entire concatenated observation table before ledger publication."}
        write_json(out / "qc.json", report)
    return report


def join_ledger(config, catalog, extracted, barcodes, output):
    """Stream a lossless, ID-checked join and validate published assignments."""
    cfg = config_read(config)
    lib = read_json(catalog)["library"]
    variants = {v["variant_id"] for v in lib["variants"]}
    observations = Path(extracted) / "observations.tsv.gz"
    assignments = Path(barcodes) / "molecule_assignments.tsv.gz"
    sources = [source_info(config, True), source_info(catalog, True), source_info(observations), source_info(assignments),
               source_info(Path(barcodes) / "mprasnakeflow_assignments.tsv.gz"),
               source_info(Path(extracted) / "qc.json", True), source_info(Path(barcodes) / "audit.json", True),
               source_info(__file__, True)]
    totals = Counter()
    with atomic_output(output) as out:
        fields = ["observed_variant_id" if f == "variant_id" else f for f in OBS_FIELDS] + [
            "corrected_barcode", "assigned_variant_id", "barcode_assignment_status"]
        with open_text(observations) as a, open_text(assignments) as b, \
                gzip.open(out / "molecule_ledger.tsv.gz", "wt", compresslevel=1, newline="") as dest:
            writer = csv.DictWriter(dest, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for obs, assigned in itertools.zip_longest(csv.DictReader(a, delimiter="\t"), csv.DictReader(b, delimiter="\t")):
                if obs is None or assigned is None or obs["read_id"] != assigned["read_id"]:
                    raise ValueError("Observation/assignment molecule counts or read IDs disagree")
                if obs["barcode"] != assigned["observed_barcode"]:
                    raise ValueError("Observed barcode changed during assignment")
                accepted = assigned["assignment_status"] in ("supported", "corrected_low_quality_substitution")
                if accepted:
                    if (assigned["variant_id"] not in variants or assigned["variant_id"] != obs["variant_id"] or
                            len(assigned["corrected_barcode"]) != lib["barcode_length"] or
                            set(assigned["corrected_barcode"]) - set("ACGT")):
                        raise ValueError("Invalid accepted barcode/genotype assignment")
                elif assigned["corrected_barcode"] or assigned["variant_id"]:
                    raise ValueError("Rejected assignment contains an accepted barcode/genotype")
                obs["observed_variant_id"] = obs.pop("variant_id")
                obs.update(corrected_barcode=assigned["corrected_barcode"], assigned_variant_id=assigned["variant_id"],
                           barcode_assignment_status=assigned["assignment_status"])
                writer.writerow(obs)
                totals["molecules"] += 1
                totals["assigned_molecules"] += accepted
                totals["analysis_eligible"] += int(obs["analysis_eligible"])
                totals["analysis_eligible_with_barcode"] += accepted and obs["analysis_eligible"] == "1"
                if totals["molecules"] % 500000 == 0:
                    print(f"Ledger joined and validated: {totals['molecules']:,} molecules", flush=True)
        qc = read_json(Path(extracted) / "qc.json")
        if totals["molecules"] != qc["counts"]["plasmid_primary"]:
            raise ValueError("Ledger is missing primary plasmid molecules")
        # Export is ordered by barcode, enabling uniqueness checking in O(1) memory.
        previous = None
        with open_text(Path(barcodes) / "mprasnakeflow_assignments.tsv.gz") as handle:
            for line in handle:
                bc, variant = line.rstrip("\n").split("\t")
                if (previous is not None and bc <= previous) or variant not in variants or len(bc) != lib["barcode_length"] or set(bc) - set("ACGT"):
                    raise ValueError("Invalid, duplicate, or unsorted assignment export")
                previous = bc
                totals["exported_barcodes"] += 1
        check_sources(sources)
        write_json(out / "validation.json", {"schema_version": "fsmpra.ledger.v1", "sources": sources,
            "counts": dict(totals), "method": cfg["barcodes"]["method"],
            "scope": qc.get("scope", "Primary plasmid molecules"),
            "observed_variant_id": "Exact catalog sequence match; confidence is recorded separately in genotype_status.",
            "assigned_variant_id": "Supported unambiguous barcode assignment consistent with the observed genotype; no rescue.",
            "limitations": ["Sequence/QC ledger only; footprint calls and occupancy statistics are not newly generated here.",
                            "Barcode coverage and effective barcode counts do not represent biological replicate counts."]})
    return dict(totals)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("extract", "join"))
    parser.add_argument("--config", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bam")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--allocated-cpus", type=int, default=4)
    parser.add_argument("--chunk-reads", type=int, default=100000)
    parser.add_argument("--extracted")
    parser.add_argument("--barcodes")
    args = parser.parse_args()
    if args.stage == "extract":
        if not args.bam:
            parser.error("--bam is required for extraction")
        extract_complete(args.config, args.catalog, args.bam, args.output, args.workers, args.allocated_cpus, args.chunk_reads)
    else:
        if not args.extracted or not args.barcodes:
            parser.error("--extracted and --barcodes are required for joining")
        join_ledger(args.config, args.catalog, args.extracted, args.barcodes, args.output)


if __name__ == "__main__":
    main()
