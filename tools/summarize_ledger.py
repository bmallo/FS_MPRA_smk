#!/usr/bin/env python3
"""Summarize a completed ledger using its small audited coverage tables."""
import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


def read_json(path):
    with open(path) as handle:
        return json.load(handle)


def summarize(run):
    run = Path(run)
    validation = read_json(run / "ledger/validation.json")
    qc = read_json(run / "extracted/qc.json")
    audit = read_json(run / "barcodes/audit.json")
    catalog = read_json(run / "catalog/catalog.json")["library"]
    variants = {v["variant_id"]: v for v in catalog["variants"]}
    groups = defaultdict(list)
    wt = {}
    with open(run / "barcodes/variant_barcode_coverage.tsv") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            kind = variants[row["variant_id"]]["class"]
            for key in ("molecules", "assigned_molecules", "unique_barcodes", "observed_barcodes"):
                row[key] = int(row[key])
            for key in ("effective_barcodes", "largest_barcode_fraction"):
                row[key] = float(row[key]) if row[key] else None
            groups[row["cohort"], kind].append(row)
            if kind == "WT":
                wt[row["cohort"]] = row
    result = []
    for (cohort, kind), rows in sorted(groups.items()):
        support = [r["unique_barcodes"] for r in rows]
        effective = [r["effective_barcodes"] for r in rows]
        result.append({"cohort": cohort, "class": kind, "designed_variants": len(rows),
            "variants_with_molecules": sum(r["molecules"] > 0 for r in rows),
            "molecules": sum(r["molecules"] for r in rows),
            "assigned_molecules": sum(r["assigned_molecules"] for r in rows),
            "sum_variant_barcode_counts": sum(support),
            "variants_with_at_least_n_supported_barcodes": {str(n): sum(x >= n for x in support) for n in (1, 5, 10, 20, 50)},
            "molecules_per_variant_min_median_max": [min(r["molecules"] for r in rows), statistics.median(r["molecules"] for r in rows), max(r["molecules"] for r in rows)],
            "supported_barcodes_per_variant_min_median_max": [min(support), statistics.median(support), max(support)],
            "effective_barcodes_per_variant_min_median_max": [min(effective), statistics.median(effective), max(effective)]})
    if sum(audit["counts"].values()) != validation["counts"]["molecules"]:
        raise ValueError("Assignment audit and validated ledger counts differ")
    if qc["counts"]["plasmid_primary"] != validation["counts"]["molecules"]:
        raise ValueError("Extraction audit and validated ledger counts differ")
    return {"schema_version": "fsmpra.ledger_summary.v1", "run": str(run.resolve()),
        "validated_counts": validation["counts"], "extraction_qc": qc["counts"],
        "molecule_assignment_status": audit["counts"], "barcode_status": audit["barcode_status_counts"],
        "assignment_evidence_molecules": audit["assignment_evidence_molecules"],
        "sequence_qualified_molecules": audit["sequence_qualified_molecules"],
        "wild_type": wt, "coverage_by_class": result,
        "interpretation": "Descriptive coverage, not power or independence. Counts across variants may share conflicting observed barcodes; only supported assignments count toward unique_barcodes. DEL/INS are intentionally excluded from analysis_eligible."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = summarize(args.run)
    with open(args.output, "x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
