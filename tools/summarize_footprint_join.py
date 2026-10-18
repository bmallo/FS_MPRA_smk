#!/usr/bin/env python3
"""Report completed footprint-join coverage from the small audited tables."""
import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


def summarize(run):
    run = Path(run)
    audit = json.loads((run / "audit.json").read_text())
    with (run / "variant_footprint_coverage.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    numeric = list(rows[0])[3:]
    for row in rows:
        for key in numeric:
            row[key] = int(row[key])
    for key in numeric:
        if sum(r[key] for r in rows) != audit["counts"].get(key, 0):
            raise ValueError("Variant coverage does not reconcile with audit: " + key)
    groups = defaultdict(list)
    for row in rows:
        groups[row["genotype_class"]].append(row)
    result = {}
    for kind, records in sorted(groups.items()):
        measured = [r["footprint_analysis_eligible"] for r in records]
        barcoded = [r["footprint_analysis_eligible_with_barcode"] for r in records]
        result[kind] = {"catalog_rows": len(records), "totals": {key: sum(r[key] for r in records) for key in numeric},
            "variants_with_measured_molecules": sum(n > 0 for n in measured),
            "variants_with_measured_barcoded_molecules": sum(n > 0 for n in barcoded),
            "measured_molecules_min_median_max": [min(measured), statistics.median(measured), max(measured)],
            "measured_barcoded_molecules_min_median_max": [min(barcoded), statistics.median(barcoded), max(barcoded)]}
    return {"schema_version": "fsmpra.footprint_join_summary.v1", "run": str(run.resolve()),
            "counts": audit["counts"], "by_genotype_class": result,
            "lowest_coverage_snvs": sorted(groups["SNV"], key=lambda r: r["footprint_analysis_eligible_with_barcode"])[:10],
            "interpretation": "Coverage only. Barcode-linked molecules are not independent barcodes or biological replicates. No occupancy effects or significance tests are computed."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = summarize(args.run)
    with open(args.output, "x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
