"""Conservative, quality-aware barcode assignment with a lossless molecule join.

This is NOT Pacybara. The exact baseline and optional one-substitution correction
provide auditable inputs/benchmarks for a future Pacybara adapter. No indel rescue
or genotype imputation is attempted. Closely related high-quality barcodes remain
distinct, even if they carry the same variant.
"""

import csv
import json
from collections import Counter
from pathlib import Path

from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .processing import SCHEMA, check_sources, config_read
from .storage import scratch_database


def assign_barcodes(config_path, catalog_path, observations_path, output):
    cfg = config_read(config_path)
    opts = cfg["barcodes"]
    method = opts["method"]
    if method not in ("exact", "quality_directional"):
        raise ValueError("Unknown barcode method")
    if (opts["min_assignment_reads"] < 1 or opts["min_seed_reads"] < 2 or opts["abundance_ratio"] < 2 or
            not 0 <= opts["max_error_quality"] < opts["min_seed_quality"] <= 93):
        raise ValueError("Invalid barcode parameters")
    lib = read_json(catalog_path)["library"]
    sources = [source_info(config_path, True), source_info(catalog_path, True), source_info(observations_path)]
    variants = {v["variant_id"]: v for v in lib["variants"]}
    counts = Counter()
    with atomic_output(output) as out, scratch_database(cfg["scratch_dir"], max_disk_gib=opts["max_disk_gib"]) as db:
        db.execute("CREATE TABLE reads(id TEXT PRIMARY KEY, barcode TEXT, variant TEXT, qualities TEXT, usable INT, eligible INT, identity_ok INT, seed_quality INT)")
        with open_text(observations_path) as handle:
            for i, row in enumerate(csv.DictReader(handle, delimiter="\t"), 1):
                usable = row["genotype_status"] == "exact_design" and row["sequence_qc"] == "1" and row["barcode_status"] == "ok"
                qs = [int(q) for q in row["barcode_qualities"].split(",") if q]
                seed_quality = usable and bool(qs) and min(qs) >= opts["min_seed_quality"]
                db.execute("INSERT INTO reads VALUES (?,?,?,?,?,?,?,?)", (row["read_id"], row["barcode"], row["variant_id"], row["barcode_qualities"], int(usable), int(row["analysis_eligible"]), int(row["genotype_status"] == "exact_design" and row["sequence_qc"] == "1"), int(seed_quality)))
                if i % 10000 == 0:
                    db.commit()
                if i % 500000 == 0:
                    print(f"Barcode evidence loaded: {i:,} molecules", flush=True)
        db.execute("CREATE INDEX reads_bc ON reads(barcode)")
        db.execute("CREATE INDEX reads_variant ON reads(variant)")
        db.execute("CREATE TABLE mapping(barcode TEXT PRIMARY KEY, corrected TEXT, variant TEXT, status TEXT, n INT, seed INT)")
        # Aggregate once; avoid a separate genotype/quality query for each of
        # millions of barcode sequences. Preserve exactly the baseline policy.
        db.execute("CREATE TABLE support AS SELECT barcode,variant,count(*) AS n,sum(seed_quality) AS seed_n FROM reads WHERE usable=1 GROUP BY barcode,variant")
        db.execute("CREATE INDEX support_bc ON support(barcode)")
        db.execute("""INSERT INTO mapping
            SELECT totals.barcode,totals.barcode,
              CASE WHEN count(s.variant)=1 THEN min(s.variant) ELSE '' END,
              CASE WHEN count(s.variant)>1 THEN 'collision'
                   WHEN sum(s.n)>=? THEN 'supported' ELSE 'insufficient_support' END,
              totals.n, CASE WHEN count(s.variant)=1 AND sum(s.n)>=? AND sum(s.seed_n)>=? THEN 1 ELSE 0 END
            FROM (SELECT barcode,count(*) AS n FROM reads GROUP BY barcode) totals
            LEFT JOIN support s ON s.barcode=totals.barcode GROUP BY totals.barcode""",
            (opts["min_assignment_reads"], opts["min_assignment_reads"], opts["min_seed_reads"]))
        db.commit()
        print("Exact barcode support aggregated", flush=True)
        write_table(out / "barcode_variant_support.tsv.gz", (
            dict(zip(["barcode", "variant_id", "supporting_molecules", "seed_quality_molecules"], row))
            for row in db.execute("SELECT * FROM support ORDER BY barcode,variant")),
            ["barcode", "variant_id", "supporting_molecules", "seed_quality_molecules"])
        corrections = 0
        if method == "quality_directional":
            # Fixed seeds, no transitive correction. Require ALL child observations
            # to agree on genotype, and all error-position qualities to be low.
            for barcode, n in db.execute("SELECT barcode,n FROM mapping WHERE status='insufficient_support'"):
                if n > opts["max_correction_reads"] or len(barcode) != lib["barcode_length"] or set(barcode) - set("ACGT"):
                    continue
                if db.execute("SELECT count(*) FROM reads WHERE barcode=? AND identity_ok=0", (barcode,)).fetchone()[0]:
                    continue
                rows = list(db.execute("SELECT variant,qualities FROM reads WHERE barcode=?", (barcode,)))
                genotypes = {v for v, _ in rows}
                if len(genotypes) != 1 or "" in genotypes:
                    continue
                qs = [[int(q) for q in qualities.split(",") if q] for _, qualities in rows]
                if any(len(q) != len(barcode) for q in qs):
                    continue
                matches = []
                for i, base in enumerate(barcode):
                    if max(q[i] for q in qs) > opts["max_error_quality"]:
                        continue
                    if any(any(q[j] < cfg["qc"]["min_base_quality"] for j in range(len(barcode)) if j != i) for q in qs):
                        continue
                    for alt in "ACGT":
                        if alt == base:
                            continue
                        parent = barcode[:i] + alt + barcode[i + 1:]
                        result = db.execute("SELECT variant,n FROM mapping WHERE barcode=? AND seed=1", (parent,)).fetchone()
                        if result and result[0] in genotypes and result[1] >= opts["abundance_ratio"] * n:
                            matches.append(parent)
                if len(matches) == 1:
                    parent = matches[0]
                    db.execute("UPDATE mapping SET corrected=?,variant=?,status='corrected_low_quality_substitution' WHERE barcode=?", (parent, next(iter(genotypes)), barcode))
                    corrections += 1
                elif len(matches) > 1:
                    db.execute("UPDATE mapping SET status='ambiguous_correction' WHERE barcode=?", (barcode,))
        db.commit()
        # Export sequences, never internal cluster IDs. Do not export error aliases
        # as extra real barcodes into RNA/DNA counting.
        with open_text(out / "mprasnakeflow_assignments.tsv.gz", "wt") as handle:
            for barcode, variant in db.execute("SELECT barcode,variant FROM mapping WHERE status='supported' ORDER BY barcode"):
                handle.write(barcode + "\t" + variant + "\n")
        fields = ["observed_barcode", "corrected_barcode", "variant_id", "status", "observed_reads", "is_seed"]
        write_table(out / "barcode_map.tsv.gz", (dict(zip(fields, r)) for r in db.execute("SELECT * FROM mapping ORDER BY barcode")), fields)
        assignment_fields = ["read_id", "observed_barcode", "corrected_barcode", "variant_id", "assignment_status"]
        with open_text(out / "molecule_assignments.tsv.gz", "wt") as handle:
            writer = csv.DictWriter(handle, fieldnames=assignment_fields, delimiter="\t")
            writer.writeheader()
            for rid, barcode, observed_variant, corrected, assigned_variant, status in db.execute(
                    "SELECT r.id,r.barcode,r.variant,m.corrected,m.variant,m.status FROM reads r JOIN mapping m ON r.barcode=m.barcode ORDER BY r.rowid"):
                if status in ("supported", "corrected_low_quality_substitution"):
                    if observed_variant != assigned_variant:
                        status = "genotype_disagreement" if observed_variant else "genotype_unresolved"
                accepted = status in ("supported", "corrected_low_quality_substitution")
                writer.writerow(dict(zip(assignment_fields, [rid, barcode, corrected if accepted else "", assigned_variant if accepted else "", status])))
                counts[status] += 1
        coverage_fields = ["variant_id", "source_variant", "cohort", "molecules", "assigned_molecules", "unique_barcodes", "effective_barcodes", "largest_barcode_fraction", "observed_barcodes"]
        coverage = []
        for variant, v in variants.items():
            for cohort, predicate in (("sequence_assigned", "1"), ("analysis_eligible", "r.eligible=1"), ("sequence_qualified", "r.identity_ok=1")):
                total = db.execute(f"SELECT count(*) FROM reads r WHERE r.variant=? AND {predicate}", (variant,)).fetchone()[0]
                groups = [n for _, n in db.execute(f"SELECT m.corrected,count(*) FROM reads r JOIN mapping m ON r.barcode=m.barcode WHERE r.variant=? AND m.variant=r.variant AND m.status IN ('supported','corrected_low_quality_substitution') AND {predicate} GROUP BY m.corrected", (variant,))]
                assigned = sum(groups)
                observed = db.execute(f"SELECT count(DISTINCT r.barcode) FROM reads r WHERE r.variant=? AND {predicate} AND length(r.barcode)=? AND r.barcode NOT GLOB '*[^ACGT]*'", (variant, lib["barcode_length"])).fetchone()[0]
                coverage.append(dict(zip(coverage_fields, [variant, v["source_variant"], cohort, total, assigned, len(groups),
                                assigned ** 2 / sum(n*n for n in groups) if groups else 0,
                                max(groups) / assigned if groups else None, observed])))
            if len(coverage) % 300 == 0:
                print(f"Variant coverage summarized: {len(coverage)//3:,} variants", flush=True)
        write_table(out / "variant_barcode_coverage.tsv", coverage, coverage_fields)
        check_sources(sources)
        report = {"schema_version": SCHEMA, "sources": sources, "method": method, "parameters": opts,
                  "counts": dict(counts), "corrected_barcode_sequences": corrections,
                  "barcode_status_counts": dict(db.execute("SELECT status,count(*) FROM mapping GROUP BY status")),
                  "sequence_qualified_molecules": db.execute("SELECT count(*) FROM reads WHERE identity_ok=1").fetchone()[0],
                  "assignment_evidence_molecules": db.execute("SELECT count(*) FROM reads WHERE usable=1").fetchone()[0],
                  "export": {"format": "headerless gzip TSV: barcode, element_id", "orientation": lib["barcode_orientation"],
                             "library_id": lib["id"], "status": "provisional; benchmark before functional production use"},
                  "limitations": ["Not Pacybara; no indel correction or genotype rescue.",
                                  "Single-read barcode evidence retained but excluded from default assignment export.",
                                  "Effective barcode count is inverse concentration, not biological replicate count.",
                                  "Barcode absence from Fiber-seq does not establish absence from the library."]}
        write_json(out / "audit.json", report)
    return report
