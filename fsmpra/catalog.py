"""Sequence-derived library catalog and allele-aware published-MPRA import."""

import csv
import hashlib
import re
from collections import Counter, defaultdict

from . import SCHEMA_VERSION, __version__
from .io import (load_manifest, new_output, read_fasta, require_columns,
                 resolve_input, source_info, write_json, write_table)


def identify_single_edit(wt, sequence, oligo_start):
    """Return canonical leftmost edit plus all sequence-equivalent edit positions.

    Positions are oligo 1-based. An insertion position is the base *after* the
    insertion; its source label contains the two flanking positions.
    """
    if sequence == wt:
        return {"class": "WT", "source_variant": "WT", "oligo_pos1": None,
                "ref": "", "alt": "", "equivalent_positions": []}
    if len(sequence) == len(wt):
        changes = [i for i, (a, b) in enumerate(zip(wt, sequence)) if a != b]
        if len(changes) != 1:
            raise ValueError("Design contains a multi-substitution sequence")
        i = changes[0]
        pos = oligo_start + i
        return {"class": "SNV", "source_variant": f"{pos}{wt[i]}>{sequence[i]}",
                "oligo_pos1": pos, "ref": wt[i], "alt": sequence[i],
                "equivalent_positions": [pos]}
    if len(sequence) == len(wt) - 1:
        candidates = [i for i in range(len(wt)) if wt[:i] + wt[i + 1:] == sequence]
        if candidates:
            pos = oligo_start + candidates[0]
            return {"class": "DEL", "source_variant": f"{pos}del{wt[candidates[0]]}",
                    "oligo_pos1": pos, "ref": wt[candidates[0]], "alt": "",
                    "equivalent_positions": [oligo_start + i for i in candidates]}
    if len(sequence) == len(wt) + 1:
        candidates = [i for i in range(len(sequence)) if sequence[:i] + sequence[i + 1:] == wt]
        if candidates:
            i = candidates[0]
            pos = oligo_start + i
            return {"class": "INS", "source_variant": f"{pos - 1}_{pos}ins{sequence[i]}",
                    "oligo_pos1": pos, "ref": "", "alt": sequence[i],
                    "equivalent_positions": [oligo_start + k for k in candidates]}
    raise ValueError("Design sequence is not WT or a single-base edit")


def source_variant_sequence(wt, label, oligo_start):
    if label == "WT":
        return wt
    match = re.fullmatch(r"(\d+)([ACGT])>([ACGT])", label)
    if match:
        pos, ref, alt = match.groups()
        i = int(pos) - oligo_start
        if not 0 <= i < len(wt) or wt[i] != ref or ref == alt:
            raise ValueError(f"Invalid reference allele: {label}")
        return wt[:i] + alt + wt[i + 1:]
    match = re.fullmatch(r"(\d+)del([ACGT])", label)
    if match:
        pos, ref = match.groups()
        i = int(pos) - oligo_start
        if not 0 <= i < len(wt) or wt[i] != ref:
            raise ValueError(f"Invalid deletion: {label}")
        return wt[:i] + wt[i + 1:]
    match = re.fullmatch(r"(\d+)_(\d+)ins([ACGT])", label)
    if match:
        left, right, alt = match.groups()
        i = int(right) - oligo_start
        if int(right) != int(left) + 1 or not 0 <= i <= len(wt):
            raise ValueError(f"Invalid insertion: {label}")
        return wt[:i] + alt + wt[i:]
    raise ValueError(f"Unsupported single-variant label: {label}")


def build_library(manifest, library):
    reference = resolve_input(manifest, library["reference_fasta"])
    seq = read_fasta(reference, library["contig"])
    digest = hashlib.sha256(seq.encode()).hexdigest()
    if digest != library["reference_sequence_sha256"]:
        raise ValueError(f"Reference checksum mismatch for {library['id']}")
    start, end = library["mutable_ref_start0"], library["mutable_ref_end0"]
    if not 0 <= start < end <= len(seq):
        raise ValueError("Invalid mutable interval")
    if library["genomic_strand"] != "+":
        raise ValueError("This adapter currently supports positive-strand libraries only")
    wt = seq[start:end]
    if set(wt) - set("ACGT"):
        raise ValueError("Mutable reference must contain only A/C/G/T")
    grouped = defaultdict(list)
    names = set()
    design_path = resolve_input(manifest, library["design_oligos"])
    with open(design_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            name, oligo = line.split()
            if name in names:
                raise ValueError(f"Duplicate oligo name {name}")
            names.add(name)
            left, right = library["design_trim_left"], library["design_trim_right"]
            sequence = oligo.upper()[left:len(oligo) - right if right else None]
            if not sequence or set(sequence) - set("ACGT"):
                raise ValueError(f"Invalid oligo {name}")
            grouped[sequence].append(name)
    variants = []
    for sequence, oligos in grouped.items():
        edit = identify_single_edit(wt, sequence, library["oligo_mutable_start1"])
        seqhash = hashlib.sha256(sequence.encode()).hexdigest()
        pos = edit["oligo_pos1"]
        offset = pos - library["oligo_mutable_start1"] if pos is not None else None
        # INS position denotes the right-flank boundary, not a substituted base.
        variants.append({"variant_id": f"{library['id']}:{seqhash}",
                         "library_id": library["id"], "sequence_sha256": seqhash,
                         "sequence": sequence, **edit,
                         "ref_pos0": start + offset if offset is not None else None,
                         "genomic_pos1": library["genomic_mutable_start1"] + offset if offset is not None else None,
                         "n_oligos": len(oligos), "oligo_names": oligos})
    variants.sort(key=lambda v: (v["oligo_pos1"] or 0, v["class"], v["source_variant"]))
    counts = dict(Counter(v["class"] for v in variants))
    if library.get("expected_unique_classes") != counts:
        raise ValueError(f"Unexpected design class counts: {counts}")
    return {**library, "reference_length": len(seq), "mutable_sequence": wt,
            "n_design_oligos": len(names), "class_counts": counts,
            "variants": variants,
            "sources": [source_info(reference, True), source_info(design_path, True)]}


def import_mpra(manifest, libraries):
    rows, audits = [], []
    for source in manifest["mpra_sources"]:
        path = resolve_input(manifest, source["path"])
        counts = Counter()
        seen = set()
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter=source["delimiter"])
            require_columns(reader, ["Chromosome", "Position", "Ref", "Alt", "Value", "P-Value", "Element", "Tags", "DNA", "RNA"], path)
            for index, row in enumerate(reader, 1):
                pos, ref, alt = int(row["Position"]), row["Ref"].upper(), row["Alt"].upper()
                chrom = "chr" + row["Chromosome"].removeprefix("chr")
                key = chrom, pos, ref, alt
                if key in seen:
                    raise ValueError(f"Duplicate MPRA allele in {path}: {key}")
                seen.add(key)
                value, p = float(row["Value"]), float(row["P-Value"])
                if not -float("inf") < value < float("inf") or not 0 <= p <= 1:
                    raise ValueError(f"Invalid MPRA value at {key}")
                matches = []
                for lib in libraries:
                    i = pos - lib["genomic_mutable_start1"]
                    if chrom != lib["genomic_chrom"] or not 0 <= i < len(lib["mutable_sequence"]):
                        continue
                    wt = lib["mutable_sequence"]
                    if wt[i] != ref:
                        raise ValueError(f"MPRA reference mismatch at {key}")
                    sequence = wt[:i] + ("" if alt == "-" else alt) + wt[i + 1:]
                    matches.extend(v["variant_id"] for v in lib["variants"] if v["sequence"] == sequence)
                if len(matches) > 1:
                    raise ValueError("MPRA matches multiple registered libraries; use an explicit join table")
                status = "matched" if matches else "outside_design"
                counts[status] += 1
                if matches and alt in "ACGT" and len(alt) == 1:
                    counts["matched_snvs"] += 1
                rows.append({"mpra_source_id": source["id"], "source_row": index,
                             "element": row["Element"], "chrom": chrom, "genomic_pos1": pos,
                             "ref": ref, "alt": alt, "effect_log2": value, "p_value": p,
                             "tags": int(row["Tags"]), "dna": float(row["DNA"]), "rna": float(row["RNA"]),
                             "variant_id": matches[0] if matches else None, "join_status": status})
        audits.append({"id": source["id"], **dict(counts), "source": source_info(path, True)})
    return rows, audits


def build_catalog(manifest_path, output):
    manifest = load_manifest(manifest_path)
    libraries = [build_library(manifest, lib) for lib in manifest["libraries"]]
    mpra, mpra_audit = import_mpra(manifest, libraries)
    samples, coverage = [], []
    for sample in manifest["samples"]:
        inputs = {key: source_info(resolve_input(manifest, value))
                  for key, value in sample["inputs"].items()}
        samples.append({**sample, "inputs": inputs})
        if "design_coverage" in sample["inputs"]:
            lib = next(lib for lib in libraries if lib["id"] == sample["library_id"])
            by_sequence = {v["sequence"]: v for v in lib["variants"]}
            path = resolve_input(manifest, sample["inputs"]["design_coverage"])
            with open(path) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                require_columns(reader, ["variant", "cls", "n_oligos", "all", "final"], path)
                seen = set()
                for row in reader:
                    sequence = source_variant_sequence(lib["mutable_sequence"], row["variant"], lib["oligo_mutable_start1"])
                    variant = by_sequence[sequence]
                    if variant["variant_id"] in seen or row["cls"] != variant["class"] or int(row["n_oligos"]) != variant["n_oligos"]:
                        raise ValueError(f"Inconsistent design coverage: {row['variant']}")
                    seen.add(variant["variant_id"])
                    for col in ("all", "pass", "strict", "all_full", "pass_full", "strict_full", "final"):
                        count = int(row[col])
                        if count < 0:
                            raise ValueError("Negative source coverage")
                        coverage.append({"sample_id": sample["id"], "variant_id": variant["variant_id"],
                                         "source_variant": row["variant"], "source_cohort": col,
                                         "count": count, "status": "historical_exploratory"})
                if len(seen) != len(lib["variants"]):
                    raise ValueError("Source coverage does not account for every design sequence")
    result = {"schema_version": SCHEMA_VERSION, "software_version": __version__,
              "manifest": source_info(manifest_path, True), "libraries": libraries,
              "samples": samples, "mpra_sources": manifest["mpra_sources"],
              "mpra_audit": mpra_audit, "mpra": mpra,
              "limitations": ["Historical coverage cohorts retain their original exploratory filters.",
                              "Published MPRA sources are kept separate; no new activity estimates or pooled p-values.",
                              "HEK293T library mapping is pending; no cross-library comparisons are made."]}
    out = new_output(output)
    write_json(out / "catalog.json", result)
    write_json(out / "manifest.snapshot.json", manifest)
    variant_fields = ["variant_id", "library_id", "source_variant", "class", "oligo_pos1", "ref_pos0", "genomic_pos1", "ref", "alt", "sequence_sha256", "sequence", "n_oligos"]
    write_table(out / "variants.tsv", ({k: v[k] for k in variant_fields} for lib in libraries for v in lib["variants"]), variant_fields)
    write_table(out / "mpra.tsv", mpra, list(mpra[0]))
    write_table(out / "historical_coverage.tsv", coverage, ["sample_id", "variant_id", "source_variant", "source_cohort", "count", "status"])
    return result
