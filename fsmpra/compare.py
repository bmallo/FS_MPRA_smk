"""Condition-aware descriptive comparisons; never pool reaction conditions."""

import csv
import itertools
import math
import statistics
from collections import Counter
from pathlib import Path

from . import SCHEMA_VERSION
from .io import load_manifest, new_output, open_text, read_json, resolve_input, source_info, write_json, write_table

COLLECTION_SCHEMA = "fsmpra.collection.v1"


def pearson(pairs):
    """Descriptive correlation, unavailable for insufficient or constant data."""
    if len(pairs) < 3:
        return None
    x, y = zip(*pairs)
    mx, my = statistics.mean(x), statistics.mean(y)
    covariance = sum((a - mx) * (b - my) for a, b in pairs)
    denominator = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return covariance / denominator if denominator else None


def processing_evidence(sample, manifest):
    """Read headers/logs, not an unbounded BAM scan or inferred wet-lab protocol."""
    import pysam
    bam_path = resolve_input(manifest, sample["inputs"]["footprints_bam"])
    with pysam.AlignmentFile(str(bam_path), "rb", threads=1) as bam:
        header = bam.header.to_dict()
    programs = sorted({(p.get("PN", p["ID"]), p.get("VN", "unknown")) for p in header.get("PG", [])
                       if p.get("PN", p["ID"]) in ("ccs", "jasmine", "fibertools-rs", "pbmm2")})
    logs = []
    for path in sorted((bam_path.parent.parent / "logs").glob("*.fiberhmm.log")):
        with path.open() as handle:
            lines = list(itertools.islice(handle, 12))
        logs.append({"source": source_info(path, True), "settings": [line.strip() for line in lines
                    if any(key in line for key in ("apply model:", "recall model:", "mode=", "min_llr=", "cores="))]})
    return {"program_versions": [{"program": p, "version": v} for p, v in programs],
            "fiberhmm_logs": logs,
            "limitation": "Logged model paths/settings do not establish historical model-file or caller-code identity."}


def check_compatible(datasets):
    first = datasets[0]
    ids, sources = set(), set()
    method_keys = ("caller", "MA_coordinates", "tracks", "min_passes", "methylation_selection",
                   "common_site_exclusion_radius_bp", "histogram_edges", "cohort")
    identity = lambda d: {(v["variant_id"], v["class"], v["oligo_pos1"], v["ref"], v["alt"]) for v in d["variants"]}
    for data in datasets:
        if data.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Unsupported measurement schema")
        if data["sample"]["id"] in ids:
            raise ValueError("Duplicate sample ID; within-sample reruns are not biological replicates")
        ids.add(data["sample"]["id"])
        bam = data["provenance"]["bam"]["path"]
        if bam in sources:
            raise ValueError("Duplicate BAM source; reruns are not independent samples")
        sources.add(bam)
        if len(identity(data)) != len(data["variants"]) or sum(v["class"] == "WT" for v in data["variants"]) != 1:
            raise ValueError("Duplicate variants or invalid WT identity")
        if any(data[k] != first[k] for k in ("library_id", "reference_sequence_sha256")):
            raise ValueError("Library/reference mismatch")
        if any(data["region"][k] != first["region"][k] for k in ("contig", "start0", "end0", "oligo_start1", "sequence")):
            raise ValueError("Analysis interval mismatch")
        if identity(data) != identity(first):
            raise ValueError("Variant/allele mismatch")
        if any(data["method"].get(k) != first["method"].get(k) for k in method_keys):
            raise ValueError("Measurement method mismatch")
        contract = data["method"].get("footprint_contract", {})
        other = first["method"].get("footprint_contract", {})
        if contract.get("format") != other.get("format") or set(contract.get("evaluated_tracks", [])) != set(other.get("evaluated_tracks", [])):
            raise ValueError("Footprint missingness contracts differ")


def sample_qc(data, min_reads):
    wt = next(v for v in data["variants"] if v["class"] == "WT")
    snvs = [v for v in data["variants"] if v["class"] == "SNV"]
    molecule_path = Path(data["provenance"]["molecules"]["path"])
    if source_info(molecule_path, True) != data["provenance"]["molecules"]:
        raise ValueError("Molecule source changed since measurement")
    classes, metrics_present = Counter(), 0
    with open_text(molecule_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["metrics_present"] == "1":
                classes[row["plasmid_class"]] += 1
                metrics_present += 1
    return {"sample_id": data["sample"]["id"], "reaction_condition_id": data["sample"].get("reaction_condition_id"),
            "sampling_fraction": data["sampling"]["fraction"], "imported_molecules": data["counts"]["imported"],
            "metrics_present": metrics_present, "plasmid_classes": dict(classes),
            "full_plasmid_fraction_of_metrics_matched_sample": classes["plasmid_full"] / metrics_present if metrics_present else None,
            "measured_molecules": data["counts"]["measured"], "WT_molecules": wt["n_reads"],
            "SNVs_total": len(snvs), "SNVs_observed": sum(v["n_reads"] > 0 for v in snvs),
            "SNVs_above_display_minimum": sum(v["n_reads"] >= min_reads for v in snvs),
            "SNV_median_molecules": statistics.median(v["n_reads"] for v in snvs) if snvs else None,
            "WT_msp_fraction_mean": wt["summary"]["msp_fraction_mean"],
            "counts": data["counts"]}


def pairwise(datasets, min_reads):
    result = []
    for a, b in itertools.combinations(datasets, 2):
        other = {v["variant_id"]: v for v in b["variants"]}
        shared = [(v, other[v["variant_id"]]) for v in a["variants"] if v["class"] == "SNV"
                  and v["n_reads"] >= min_reads and other[v["variant_id"]]["n_reads"] >= min_reads]
        ca, cb = a["sample"].get("reaction_condition_id"), b["sample"].get("reaction_condition_id")
        relation = "same_stated_reaction" if ca and ca == cb else "different_reaction_conditions" if ca and cb else "reaction_unknown"
        entry = {"sample_a": a["sample"]["id"], "sample_b": b["sample"]["id"],
                 "comparison_type": relation, "SNVs_above_minimum_in_both": len(shared)}
        for metric in ("msp_common_site_delta", "tf_small_common_site_delta"):
            pairs = [(x["summary"][metric], y["summary"][metric]) for x, y in shared
                     if x["summary"][metric] is not None and y["summary"][metric] is not None]
            entry[metric + "_n"] = len(pairs)
            entry[metric + "_pearson"] = pearson(pairs)
        result.append(entry)
    return result


def recovery_evidence(data):
    """Optional full-import QC; older pilot bundles remain readable."""
    path = Path(data["provenance"]["molecules"]["path"]).parent / "recovery.json"
    if not path.exists():
        return None
    recovery = read_json(path)
    if (recovery["sample_id"] != data["sample"]["id"] or recovery["fraction"] != data["sampling"]["fraction"]
            or sum(g["molecules"] for g in recovery["by_plasmid_class"].values()) != data["counts"]["imported"]):
        raise ValueError("Recovery summary does not match measurement population")
    return {"source": source_info(path, True), "summary": recovery}


def compare_pilots(paths, manifest_path, output, min_reads=10, include_processing=True):
    if min_reads < 1 or len(paths) < 2:
        raise ValueError("Comparison needs at least two samples and a positive display minimum")
    datasets = [read_json(p) for p in paths]
    check_compatible(datasets)
    manifest = load_manifest(manifest_path)
    samples = {s["id"]: s for s in manifest["samples"]}
    libraries = {l["id"]: l for l in manifest["libraries"]}
    processing, rows, conditions = {}, [], {}
    for data in datasets:
        sample = samples[data["sample"]["id"]]
        lib = libraries[sample["library_id"]]
        condition = sample.get("reaction_condition_id")
        reaction = sample.get("reaction", {})
        parameters = tuple(reaction.get(k) for k in ("duration_minutes", "hia5_volume_uL", "cells"))
        if condition and condition in conditions and conditions[condition] != parameters:
            raise ValueError("One reaction condition ID has conflicting experimental parameters")
        if condition:
            conditions[condition] = parameters
        if data["library_id"] != lib["id"] or data["reference_sequence_sha256"] != lib["reference_sequence_sha256"]:
            raise ValueError("Metadata manifest library/reference mismatch")
        if (data["region"]["start0"], data["region"]["end0"], data["region"]["oligo_start1"]) != (lib["mutable_ref_start0"], lib["mutable_ref_end0"], lib["oligo_mutable_start1"]):
            raise ValueError("Metadata manifest interval mismatch")
        if source_info(resolve_input(manifest, sample["inputs"]["footprints_bam"])) != data["provenance"]["bam"]:
            raise ValueError("Metadata manifest BAM/source mismatch")
        # Amend experimental annotations in the new collection only. Never
        # rewrite the source measurements or their historical method settings.
        data["sample"].update({k: sample[k] for k in ("reaction_condition_id", "reaction", "protocol_equivalence") if k in sample})
        data["region"]["status"] = lib["analysis_interval_status"]
        data["provenance"]["metadata_manifest"] = source_info(manifest_path, True)
        data["limitations"] = [s for s in data["limitations"] if "promoter boundary selection remains provisional" not in s]
        data["limitations"].append("Reaction conditions differ across this collection; runs are not pooled or treated as exchangeable replicates.")
        if include_processing:
            processing[sample["id"]] = processing_evidence(sample, manifest)
        for v in data["variants"]:
            rows.append({"sample_id": sample["id"], "reaction_condition_id": sample.get("reaction_condition_id"),
                         "variant_id": v["variant_id"], "source_variant": v["source_variant"], "class": v["class"],
                         "oligo_pos1": v["oligo_pos1"], "n_reads": v["n_reads"],
                         "msp_common_site_delta": v["summary"]["msp_common_site_delta"],
                         "tf_small_common_site_delta": v["summary"]["tf_small_common_site_delta"]})
    report = {"display_minimum_reads": min_reads, "sample_qc": [sample_qc(d, min_reads) for d in datasets],
              "pairwise": pairwise(datasets, min_reads), "processing_evidence": processing,
              "recovery": {d["sample"]["id"]: recovery_evidence(d) for d in datasets},
              "limitations": ["Descriptive comparison; no pooling, hypothesis tests, or biological-replicate uncertainty estimates.",
                              "Reaction assignments come from the supplied metadata manifest. Unrecorded preparation details remain unknown.",
                              "Reaction, run, software, and biological effects cannot be separated where a reaction condition has only one biological preparation.",
                              "Correlations use the same SNVs passing the display minimum in both samples; small or low-depth intersections can be unstable.",
                              "Identical logged FiberHMM settings do not establish identical upstream modification calls or historical model files.",
                              "Coverage and full-plasmid fractions describe the sampled genotype-ledger population, not all sequenced or transfected molecules."],
              "provenance": {"measurements": [source_info(p, True) for p in paths], "metadata_manifest": source_info(manifest_path, True)}}
    collection = {"schema_version": COLLECTION_SCHEMA, "datasets": datasets, "comparison": report}
    out = new_output(output)
    write_json(out / "collection.json", collection)
    write_json(out / "comparison.json", report)
    write_table(out / "allele_comparison.tsv", rows, list(rows[0]))
    flat = [{k: v for k, v in row.items() if k not in ("plasmid_classes", "counts")} for row in report["sample_qc"]]
    write_table(out / "sample_qc.tsv", flat, list(flat[0]))
    return {"samples": len(datasets), "sample_qc": flat, "pairwise": report["pairwise"], "output": str(out)}
