"""Bounded barcode benchmarks; never changes production assignment defaults.

Pacybara is executed from a pinned external checkout, not reimplemented here.
Its genotype input is restricted to confidently observed WT/single substitutions.
"""

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
import random
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

from .io import atomic_output, open_text, read_json, read_fasta, source_info, write_json, write_table
from .processing import OBS_FIELDS, check_sources, config_read, extract, extract_interval, prepare
from .barcodes import assign_barcodes

PACYBARA_REVISION = "b5e39bd0e2731283eee2ee88c0e0852ddf83fe33"


def hash_fraction(name, seed):
    return int.from_bytes(hashlib.sha256((seed + "\0" + name).encode()).digest()[:8], "big") / 2**64


def sample_bam(config, bam_path, output, fraction=0.01, seed="barcode-benchmark-v1"):
    """One indexed contig traversal; stable independent inclusion by molecule ID."""
    import pysam
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    cfg = config_read(config)
    sources = [source_info(config, True), source_info(bam_path), source_info(str(bam_path) + ".bai")]
    counts = Counter()
    start = time.monotonic()
    with atomic_output(output) as out:
        prepare(config, out / "catalog")
        with pysam.AlignmentFile(bam_path, "rb", threads=1) as src, \
                pysam.AlignmentFile(out / "sample.bam", "wb", template=src, threads=1) as dest:
            for read in src.fetch(cfg["library"]["contig"]):
                counts["target_alignment_records"] += 1
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                counts["primary_plasmids"] += 1
                if hash_fraction(read.query_name, seed) < fraction:
                    dest.write(read)
                    counts["sampled_primary_plasmids"] += 1
                if counts["primary_plasmids"] % 1000000 == 0:
                    print(dict(counts), flush=True)
        print("Sampling complete; extracting sequences and QC: " + json.dumps(dict(counts)), flush=True)
        extract(config, out / "catalog/catalog.json", out / "sample.bam", out / "extracted")
        # extract provenance contains the temporary publication path; record the
        # stable final location alongside the original immutable input sources.
        check_sources(sources)
        write_json(out / "sampling.json", {"sources": sources, "seed": seed,
            "fraction": fraction, "counts": dict(counts), "seconds": time.monotonic() - start,
            "selection": "SHA256(seed + NUL + read_id), first 64 bits / 2**64 < fraction",
            "final_sample_bam": str(Path(output).resolve() / "sample.bam"),
            "limitation": "Read subsampling reduces family depth; rare-family conclusions need full-depth follow-up."})


def one_edit_neighborhood(barcode):
    result = {barcode}
    for i, base in enumerate(barcode):
        result.add(barcode[:i] + barcode[i+1:])
        for alt in "ACGT":
            if alt != base:
                result.add(barcode[:i] + alt + barcode[i+1:])
    for i in range(len(barcode) + 1):
        for alt in "ACGT":
            result.add(barcode[:i] + alt + barcode[i:])
    return result


def family_bam(config, bam_path, discovery, output, families_per_stratum=200):
    """Collect full-depth edit-one neighborhoods; selection is discovery-informed."""
    import pysam
    if not 1 <= families_per_stratum <= 1000:
        raise ValueError("Choose 1–1000 families per discovery stratum")
    cfg = config_read(config)
    lib = cfg["library"]
    discovery = Path(discovery)
    obs = discovery / "extracted/observations.tsv.gz"
    sources = [source_info(config, True), source_info(bam_path), source_info(str(bam_path) + ".bai"), source_info(obs, True)]
    groups = defaultdict(list)
    for r in load_rows(obs):
        if (r["genotype_class"] in ("WT", "SNV") and r["genotype_status"] == "exact_design" and
                r["sequence_qc"] == "1" and r["barcode_status"] == "ok" and
                min(map(int, r["barcode_qualities"].split(","))) >= 30):
            groups[r["barcode"]].append(r)
    seeds = []
    for stratum in ("repeated", "singleton"):
        eligible = [bc for bc, rows in groups.items() if len({r["variant_id"] for r in rows}) == 1 and
                    (len(rows) >= 2 if stratum == "repeated" else len(rows) == 1)]
        for bc in sorted(eligible, key=lambda b: hash_fraction(b, "barcode-family-v1"))[:families_per_stratum]:
            seeds.append({"barcode": bc, "discovery_stratum": stratum, "discovery_q30_reads": len(groups[bc]),
                          "discovery_variant_id": groups[bc][0]["variant_id"]})
    if not seeds:
        raise ValueError("No eligible discovery barcodes")
    neighborhood = set().union(*(one_edit_neighborhood(s["barcode"]) for s in seeds))
    counts = Counter()
    started = time.monotonic()
    with atomic_output(output) as out:
        prepare(config, out / "catalog")
        reference = read_fasta(out / "catalog/plasmid.fa", lib["contig"])
        write_table(out / "seeds.tsv", seeds, ["barcode", "discovery_stratum", "discovery_q30_reads", "discovery_variant_id"])
        with pysam.AlignmentFile(bam_path, "rb", threads=1) as src, \
                pysam.AlignmentFile(out / "sample.bam", "wb", template=src, threads=1) as dest:
            for read in src.fetch(lib["contig"]):
                counts["target_alignment_records"] += 1
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                counts["primary_plasmids"] += 1
                bc, _, status = extract_interval(read, lib["barcode_start0"], lib["barcode_end0"], reference)
                if lib["barcode_orientation"] == "reverse_complement":
                    bc = bc.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
                if status == "ok" and bc in neighborhood:
                    dest.write(read)
                    counts["sampled_primary_plasmids"] += 1
                if counts["primary_plasmids"] % 1000000 == 0:
                    print(dict(counts), flush=True)
        print("Family collection complete; extracting sequences and QC: " + json.dumps(dict(counts)), flush=True)
        extract(config, out / "catalog/catalog.json", out / "sample.bam", out / "extracted")
        check_sources(sources)
        write_json(out / "sampling.json", {"sources": sources, "counts": dict(counts),
            "selection": "Complete anchored edit-distance <=1 barcode neighborhoods around deterministic discovery seeds",
            "families_per_stratum_requested": families_per_stratum, "seeds": len(seeds),
            "neighborhood_sequences": len(neighborhood), "seconds": time.monotonic() - started,
            "limitation": "Discovery-informed targeted families, not population-representative; neighborhoods exclude >1-edit errors relative to every seed and unreadable flanks."})


def edit_distance(a, b):
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(current[-1] + 1, previous[j] + 1, previous[j-1] + (x != y)))
        previous = current
    return previous[-1]


def nearby_pairs(barcodes, radius=2):
    """Complete Levenshtein candidates via deletion signatures, then exact DP.

    Radius two supplies Pacybara's complete-link safeguard for maxDiff=1.
    No quadratic all-barcode comparison and no external aligner heuristic.
    """
    if radius not in (1, 2):
        raise ValueError("This bounded adapter supports radius 1 or 2 only")
    index = defaultdict(list)
    for bc in sorted(set(barcodes)):
        signatures = {bc}
        for n in range(1, radius + 1):
            for positions in itertools.combinations(range(len(bc)), n):
                signatures.add("".join(x for i, x in enumerate(bc) if i not in positions))
        candidates = set()
        for signature in signatures:
            candidates.update(index.get(signature, ()))
        for other in sorted(candidates):
            if abs(len(other) - len(bc)) <= radius:
                distance = edit_distance(bc, other)
                if distance <= radius:
                    yield other, bc, distance
        for signature in signatures:
            index[signature].append(bc)


def genotype_tokens(lib):
    wt = lib["mutable_sequence"]
    tokens = {}
    for v in lib["variants"]:
        if v["class"] == "WT":
            tokens[v["variant_id"]] = "="
        elif v["class"] == "SNV":
            diffs = [(i, a, b) for i, (a, b) in enumerate(zip(wt, v["sequence"]), 1) if a != b]
            if len(diffs) != 1 or len(wt) != len(v["sequence"]):
                raise ValueError("Invalid catalog SNV")
            i, a, b = diffs[0]
            tokens[v["variant_id"]] = f"{a}{i}{b}"
    return tokens


def run_pacybara(rows, lib, checkout, output, max_pairs=3000000):
    """Run unmodified pinned core; adapt inputs and reject ambiguous exports.

    Genotype Q20 is a conservative lower bound from the common cohort filter.
    minQual=19 retains these trusted genotypes (the upstream comparison is >).
    Barcode quality is NOT used by upstream clustering. Synthetic FASTQ quality
    lines are deliberately placeholders because that reader ignores them.
    """
    checkout, output = Path(checkout).resolve(), Path(output)
    revision = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip()
    if revision != PACYBARA_REVISION or subprocess.check_output(
            ["git", "-C", str(checkout), "status", "--porcelain"], text=True).strip():
        raise ValueError("Pacybara checkout must be clean and match the pinned revision")
    output.mkdir()
    tokens = genotype_tokens(lib)
    reverse = {v: k for k, v in tokens.items()}
    groups = defaultdict(list)
    aliases = {}
    for i, row in enumerate(rows):
        alias = f"r{i:09d}"
        aliases[alias] = row
        groups[row["barcode"]].append(alias)
    write_table(output / "read_ids.tsv.gz", ({"pacybara_id": k, "read_id": r["read_id"]} for k, r in aliases.items()), ["pacybara_id", "read_id"])
    with gzip.open(output / "preclusters.fastq.gz", "wt") as handle:
        for bc, ids in sorted(groups.items()):
            handle.write(f"@{'|'.join(ids)} size={len(ids)}\n{bc}\n+\n{'I' * len(bc)}\n")
    with gzip.open(output / "genotypes.csv.gz", "wt") as handle:
        for alias, row in aliases.items():
            token = tokens[row["variant_id"]]
            handle.write(alias + "," + ("=" if token == "=" else token + ":20") + "\n")
    pairs = sum(len(ids) * (len(ids)-1) // 2 for ids in groups.values())
    if pairs > max_pairs:
        raise ValueError("Pacybara exact-family pair budget exceeded; reduce sample size")
    edge_count = 0
    with gzip.open(output / "distances.csv.gz", "wt") as handle:
        handle.write("read1,read2,barcode1,barcode2,distance\n")
        for a, b, distance in nearby_pairs(groups):
            pairs += len(groups[a]) * len(groups[b])
            if pairs > max_pairs:
                raise ValueError("Pacybara expanded read-pair budget exceeded; reduce sample size")
            edge_count += 1
            handle.write(f"{'|'.join(groups[a])},{'|'.join(groups[b])},{a},{b},{distance}\n")
    argv = [sys.executable, str(checkout / "src/pacybara_cluster.py"),
            "--maxDiff", "1", "--minQual", "19", "--minMatches", "1", "--minJaccard", "0.2",
            "--out", str(output / "clusters.csv.gz"), str(output / "distances.csv.gz"),
            str(output / "genotypes.csv.gz"), str(output / "preclusters.fastq.gz")]
    # Upstream can invoke MUSCLE for a tied >2-read consensus. Do not substitute
    # a different consensus silently: any missing dependency fails this run.
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PYTHONHASHSEED="0")
    start = time.monotonic()
    with open(output / "core.log", "w") as handle:
        subprocess.run(argv, check=True, stdout=handle, stderr=handle, env=env)
    mappings = defaultdict(set)
    assignments = {}
    seen = set()
    clusters = []
    with gzip.open(output / "clusters.csv.gz", "rt") as handle:
        for cluster in csv.DictReader(handle):
            ids = cluster["reads"].split("|")
            if int(cluster["size"]) != len(ids) or seen.intersection(ids):
                raise ValueError("Duplicate IDs or inconsistent cluster size")
            seen.update(ids)
            variant = reverse.get(cluster["geno"], "")
            # Preserve raw clusters even when final library export is ambiguous.
            bc = cluster["virtualBarcode"]
            clusters.append((bc, variant, ids))
            mappings[bc].add(variant)
    if seen != set(aliases):
        raise ValueError("Pacybara lost or introduced molecules")
    export = {}
    for bc, variant, ids in clusters:
        accepted = (len(ids) >= 2 and len(bc) == lib["barcode_length"] and
                    not set(bc) - set("ACGT") and variant and mappings[bc] == {variant})
        if accepted:
            export[bc] = variant
        for alias in ids:
            row = aliases[alias]
            # Same conservative no-genotype-rescue rule as the baseline.
            assignments[row["read_id"]] = (bc, variant) if accepted and row["variant_id"] == variant else ("", "")
    write_table(output / "molecule_assignments.tsv.gz", (
        {"read_id": rid, "corrected_barcode": bc, "variant_id": variant}
        for rid, (bc, variant) in sorted(assignments.items())), ["read_id", "corrected_barcode", "variant_id"])
    with gzip.open(output / "mprasnakeflow_assignments.tsv.gz", "wt") as handle:
        for bc, variant in sorted(export.items()):
            handle.write(f"{bc}\t{variant}\n")
    write_json(output / "adapter.json", {"revision": revision,
        "core_source": source_info(checkout / "src/pacybara_cluster.py", True),
        "argv": argv, "barcode_edges": edge_count, "expanded_pairs": pairs,
        "seconds_core_only": time.monotonic() - start, "reads": len(rows),
        "raw_clusters": len(clusters), "exported_barcodes": len(export),
        "ambiguous_consensus_barcodes": sum(len(v) > 1 for v in mappings.values()),
        "scope": "Unmodified clustering core with shared LDLR extraction and exhaustive edit-distance <=2 adapter; not full Pacybara workflow.",
        "parameters": {"maxDiff": 1, "minQual": 19, "minMatches": 1, "minJaccard": 0.2},
        "genotype_quality": "Q20 lower bound from trusted WT/SNV cohort; barcode qualities unused by upstream core",
        "export_policy": "N15 ACGT, >=2 reads, unambiguous consensus genotype, reject read-genotype disagreement"})
    return assignments, export


def load_rows(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def benchmark_methods(rows, config, catalog, checkout, output):
    output.mkdir()
    cfg = config_read(config)
    lib = read_json(catalog)["library"]
    write_table(output / "observations.tsv.gz", rows, OBS_FIELDS)
    methods = {}
    for method in ("exact", "quality_directional"):
        print(f"{output.name}: {method}, {len(rows)} molecules", flush=True)
        cfg["barcodes"]["method"] = method
        path = output / (method + ".json")
        write_json(path, cfg)
        assign_barcodes(path, catalog, output / "observations.tsv.gz", output / method)
        assignments = {r["read_id"]: (r["corrected_barcode"], r["variant_id"])
                       for r in load_rows(output / method / "molecule_assignments.tsv.gz")}
        with gzip.open(output / method / "mprasnakeflow_assignments.tsv.gz", "rt") as handle:
            export = dict(line.rstrip().split("\t") for line in handle)
        methods[method] = (assignments, export)
    print(f"{output.name}: Pacybara core, {len(rows)} molecules", flush=True)
    methods["pacybara_core_ed1"] = run_pacybara(rows, lib, checkout, output / "pacybara_core_ed1")
    coverage, merges = [], []
    support = Counter((r["barcode"], r["variant_id"]) for r in rows
                      if r["barcode_status"] == "ok" and min(map(int, r["barcode_qualities"].split(","))) >= 30)
    for method, (assignments, _) in methods.items():
        variant_counts = defaultdict(Counter)
        transitions = Counter()
        for r in rows:
            bc, variant = assignments[r["read_id"]]
            if bc:
                variant_counts[variant][bc] += 1
                if bc != r["barcode"]:
                    transitions[r["barcode"], bc, variant] += 1
        for variant, counts in sorted(variant_counts.items()):
            n = sum(counts.values())
            coverage.append({"method": method, "variant_id": variant, "assigned_molecules": n,
                "unique_barcodes": len(counts), "effective_barcodes": n*n/sum(x*x for x in counts.values()),
                "largest_barcode_fraction": max(counts.values())/n})
        for (observed, corrected, variant), n in sorted(transitions.items()):
            merges.append({"method": method, "observed_barcode": observed, "corrected_barcode": corrected,
                "variant_id": variant, "changed_molecules": n, "observed_q30_reads": support[observed, variant],
                "parent_q30_reads": support[corrected, variant]})
    write_table(output / "variant_coverage.tsv", coverage,
                ["method", "variant_id", "assigned_molecules", "unique_barcodes", "effective_barcodes", "largest_barcode_fraction"])
    write_table(output / "corrections.tsv", merges,
                ["method", "observed_barcode", "corrected_barcode", "variant_id", "changed_molecules", "observed_q30_reads", "parent_q30_reads"])
    return methods


def synthetic_cases(lib, repeats=20):
    """Adversarial known-truth cases, not an empirical error-rate simulation."""
    rng = random.Random(20261016)
    wt = next(v for v in lib["variants"] if v["class"] == "WT")
    snv = next(v for v in lib["variants"] if v["class"] == "SNV")
    rows, truth = [], {}

    def add(bc, variant, n, family, case, quality=None):
        for _ in range(n):
            rid = f"synthetic_{len(rows):06d}"
            qs = quality or [40] * len(bc)
            row = dict.fromkeys(OBS_FIELDS, "")
            row.update(read_id=rid, sample_id="synthetic", library_id=lib["id"],
                observed_sequence=variant["sequence"], variant_id=variant["variant_id"],
                genotype_class=variant["class"], genotype_status="exact_design", barcode=bc,
                barcode_qualities=",".join(map(str, qs)),
                barcode_status="length_mismatch" if len(bc) != 15 else "low_base_quality" if min(qs) < 20 else "ok",
                sequence_qc="1", analysis_eligible="1")
            rows.append(row)
            truth[rid] = {"barcode": family, "variant_id": variant["variant_id"], "case": case}

    for case in ("low_quality_substitution", "high_quality_substitution", "insertion", "deletion",
                 "genuine_neighbor_same_variant", "genuine_neighbor_other_variant", "barcode_collision", "ambiguous_parent"):
        for rep in range(repeats):
            bc = "".join(rng.choice("ACGT") for _ in range(15))
            alt = next(x for x in "ACGT" if x != bc[7])
            child = bc[:7] + alt + bc[8:]
            variant = wt if rep % 2 == 0 else snv
            other = snv if variant == wt else wt
            add(bc, variant, 15, bc, case)
            if case == "barcode_collision":
                add(bc, other, 3, bc, case)
            elif case == "ambiguous_parent":
                alt2 = next(x for x in "ACGT" if x not in (bc[7], alt))
                parent2 = bc[:7] + alt2 + bc[8:]
                add(parent2, variant, 15, parent2, case)
                add(child, variant, 1, bc, case, [40]*7 + [10] + [40]*7)
            elif case.startswith("genuine_neighbor"):
                add(child, other if case.endswith("other_variant") else variant, 3, child, case)
            else:
                if case == "insertion":
                    child = bc[:7] + alt + bc[7:]
                elif case == "deletion":
                    child = bc[:7] + bc[8:]
                q = [40] * len(child)
                if case == "low_quality_substitution":
                    q[7] = 10
                add(child, variant, 1, bc, case, q)
    return rows, truth


def summarize(rows, methods, truth=None, heldout=None):
    summaries = []
    baseline = methods.get("exact", ({}, {}))[0]
    for method, (assigned, export) in methods.items():
        groups = defaultdict(set)
        counts = Counter()
        anchors = Counter((r["barcode"], r["variant_id"]) for r in rows if
                          r["barcode_status"] == "ok" and min(map(int, r["barcode_qualities"].split(","))) >= 30)
        anchor_merges = set()
        by_case = defaultdict(Counter)
        for row in rows:
            rid = row["read_id"]
            bc, variant = assigned[rid]
            counts["input_molecules"] += 1
            if bc:
                counts["assigned_molecules"] += 1
                counts["changed_barcode_molecules"] += bc != row["barcode"]
                counts["variant_disagreements"] += variant != row["variant_id"]
                if baseline:
                    was_assigned = bool(baseline[rid][0])
                    counts["newly_assigned_vs_exact"] += not was_assigned
                    counts["changed_previously_assigned_molecules"] += was_assigned and bc != baseline[rid][0]
                if bc != row["barcode"] and anchors[row["barcode"], variant] >= 3 and anchors[bc, variant] >= 3:
                    anchor_merges.add((row["barcode"], bc, variant))
            if truth:
                t = truth[rid]
                case = by_case[t["case"]]
                case["input_molecules"] += 1
                if bc:
                    case["assigned_molecules"] += 1
                    correct = bc == t["barcode"] and variant == t["variant_id"]
                    case["correct_assignments"] += correct
                    case["incorrect_assignments"] += not correct
                    counts["incorrect_assignments"] += not correct
                    groups[bc].add((t["barcode"], t["variant_id"]))
        counts["exported_barcodes"] = len(export)
        counts["merges_between_q30_supported_barcodes"] = len(anchor_merges)
        if truth:
            counts["merged_truth_families_clusters"] = sum(len(v) > 1 for v in groups.values())
        if heldout is not None:
            counts["heldout_molecules"] = len(heldout)
            for row in heldout:
                # Freeze export; NEVER use held-out genotype to select a match.
                predicted = export.get(row["barcode"])
                if predicted:
                    counts["heldout_exact_export_matches"] += 1
                    counts["heldout_variant_disagreements"] += predicted != row["variant_id"]
        summaries.append({"method": method, **dict(counts), "synthetic_cases": dict(by_case)})
    return summaries


def split_holdout(rows, discovery_ids=()):
    discovery_ids = set(discovery_ids)
    train, heldout = [], []
    for row in rows:
        (heldout if row["read_id"] not in discovery_ids and
         hash_fraction(row["read_id"], "barcode-holdout-v1") < 0.2 else train).append(row)
    return train, heldout


def run_benchmark(config, sample, checkout, output):
    sample, output = Path(sample), Path(output)
    catalog = sample / "catalog/catalog.json"
    lib = read_json(catalog)["library"]
    cfg = config_read(config)
    extraction_qc = read_json(sample / "extracted/qc.json")
    if extraction_qc["qc"]["min_base_quality"] < 20:
        raise ValueError("Adapter genotype lower bound requires extraction min_base_quality >=20")
    if lib["id"] != cfg["library"]["id"] or lib["barcode_length"] != 15:
        raise ValueError("This benchmark requires the matching library and N15 barcodes")
    sources = [source_info(config, True), source_info(catalog, True), source_info(sample / "extracted/observations.tsv.gz", True),
               source_info(sample / "extracted/qc.json", True),
               source_info(__file__, True), source_info(Path(__file__).with_name("barcodes.py"), True),
               source_info(Path(__file__).with_name("processing.py"), True)]
    sampling = read_json(sample / "sampling.json")
    discovery_ids = set()
    if "seeds" in sampling:
        # Family selection consulted discovery genotypes. Keep ALL discovery
        # molecules out of validation, even when their hash would select them.
        discovery_source = sampling["sources"][-1]
        check_sources([discovery_source])
        sources.append(discovery_source)
        discovery_ids = {r["read_id"] for r in load_rows(discovery_source["path"])}
    all_rows = load_rows(sample / "extracted/observations.tsv.gz")
    rows = [r for r in all_rows if r["genotype_class"] in ("WT", "SNV") and
            r["genotype_status"] == "exact_design" and r["sequence_qc"] == "1" and
            r["barcode_status"] in ("ok", "low_base_quality", "length_mismatch") and
            14 <= len(r["barcode"]) <= 16 and not set(r["barcode"]) - set("ACGT") and
            len(r["barcode_qualities"].split(",")) == len(r["barcode"])]
    rows.sort(key=lambda r: r["read_id"])
    if len({r["read_id"] for r in all_rows}) != len(all_rows):
        raise ValueError("Duplicate input molecule ID")
    with atomic_output(output) as out:
        write_json(out / "config.json", config_read(config))
        synthetic, truth = synthetic_cases(lib)
        write_json(out / "synthetic_truth.json", truth)
        synthetic_methods = benchmark_methods(synthetic, config, catalog, checkout, out / "synthetic")
        print("Synthetic benchmark complete", flush=True)
        real_methods = benchmark_methods(rows, config, catalog, checkout, out / "real")
        print("Real-data benchmark complete", flush=True)
        train, heldout = split_holdout(rows, discovery_ids)
        write_table(out / "heldout_observations.tsv.gz", heldout, OBS_FIELDS)
        train_methods = benchmark_methods(train, config, catalog, checkout, out / "training")
        write_table(out / "heldout_predictions.tsv.gz", (
            {"method": method, "read_id": r["read_id"], "barcode": r["barcode"],
             "observed_variant_id": r["variant_id"], "predicted_variant_id": export.get(r["barcode"], ""),
             "outcome": "unmatched" if r["barcode"] not in export else
                        "agree" if export[r["barcode"]] == r["variant_id"] else "disagree"}
            for method, (_, export) in train_methods.items() for r in heldout),
            ["method", "read_id", "barcode", "observed_variant_id", "predicted_variant_id", "outcome"])
        report = {"schema_version": "fsmpra.barcode_benchmark.v1", "sources": sources,
            "sampling": sampling, "sampled_primary_plasmids": len(all_rows),
            "discovery_molecules_forced_into_training": sum(r["read_id"] in discovery_ids for r in rows),
            "comparison_cohort_molecules": len(rows), "comparison_unique_observed_barcodes": len({r['barcode'] for r in rows}),
            "comparison_cohort": "High-confidence WT/SNV genotypes and sequence QC, anchored ACGT barcodes length 14–16; no chromatin-outcome filtering",
            "all_sample_barcode_status": dict(Counter(r["barcode_status"] for r in all_rows)),
            "cohort_barcode_status": dict(Counter(r["barcode_status"] for r in rows)),
            "real": summarize(rows, real_methods), "synthetic": summarize(synthetic, synthetic_methods, truth),
            "holdout": summarize(train, train_methods, heldout=heldout),
            "limitations": ["Read-sampling reduces family depth; targeted-family sampling enriches discovery families. Neither estimates full-library correction yield.",
                "Real barcode truth is unknown; assignment yield is not accuracy. Training agreement is gated, not validation.",
                "Held-out genotypes never enter training. Exact matching against frozen representative export tests link consistency, not error correction or biological independence.",
                "Synthetic cases are adversarial known-truth tests, not estimates of real error prevalence.",
                "Pacybara core only, maxDiff=1 with exact radius-two distances; indel genotypes and maxDiff=2 not evaluated.",
                "MPRAsnakeflow export barcode orientation still needs validation against future functional reads."]}
        check_sources(sources)
        write_json(out / "report.json", report)
        print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)
    p = sub.add_parser("sample")
    p.add_argument("--config", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--fraction", type=float, default=0.01)
    p.add_argument("--seed", default="barcode-benchmark-v1")
    p = sub.add_parser("families")
    p.add_argument("--config", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--discovery", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--families-per-stratum", type=int, default=200)
    p = sub.add_parser("run")
    p.add_argument("--config", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--checkout", required=True)
    p.add_argument("--output", required=True)
    args = vars(parser.parse_args())
    stage = args.pop("stage")
    if stage in ("sample", "families"):
        args["bam_path"] = args.pop("bam")
        (sample_bam if stage == "sample" else family_bam)(**args)
    else:
        run_benchmark(**args)


if __name__ == "__main__":
    main()
