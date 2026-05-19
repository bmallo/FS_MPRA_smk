#!/usr/bin/env python3
"""
h5_to_json.py — Convert a Fiber-seq MPRA pipeline {sample}.h5 into a JSON
file consumable by the dark-mode dashboard (data.js shape).

Usage:
    python h5_to_json.py path/to/sample.h5 > data.json
    # then in index.html replace the inline data.js with:
    # <script>fetch('data.json').then(r=>r.json()).then(d=>{window.MPRA={...d}})</script>

Output JSON shape mirrors window.MPRA in data.js:
    {
      "samples": { "<name>": <Sample> },
      "sampleOrder": ["<name>", ...]
    }

Each Sample exposes:
    metadata, ref_sequence, wt_occupancy, ground_truth_nc,
    variants[], motifs[], density_tracks, cooccupancy[],
    tf_candidates[], mpra[]

Dependencies: h5py, numpy. (No pandas required.)
"""

from __future__ import annotations
import sys
import json
import argparse
import h5py
import numpy as np
from pathlib import Path


def _to_pylist(v):
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return v.item()
    if isinstance(v, bytes):
        return v.decode("utf-8", errors="replace")
    return v


def _attrs(group) -> dict:
    return {k: _to_pylist(group.attrs[k]) for k in group.attrs.keys()}


def _round_list(lst, precision):
    return [round(x, precision) for x in lst]


def convert(h5_path: Path, sample_name: str | None = None,
            mpra_tsv: Path | None = None, precision: int = 4) -> dict:
    with h5py.File(h5_path, "r") as f:
        meta = _attrs(f["metadata"])
        if sample_name is None:
            sample_name = meta.get("sample_name") or h5_path.stem
        meta["sample_name"] = sample_name
        bins = meta["bin_labels"]

        # WT occupancy per bin.
        # Pipeline only writes this under shared-null calibration; for
        # per-variant nulls the group is absent. Reconstruct post-hoc as
        # mean(variant_occ - delta_obs) across variants (accumulated below).
        wt = {}
        _wt_acc = {b: None for b in bins}   # filled during variant loop
        if "wt_occupancy" in f:
            wt = {b: f[f"wt_occupancy/{b}"][:].tolist() for b in bins}

        # Per-variant summary arrays indexed by variant_idx (pipeline v0.2+)
        _summary = {}
        if "summary" in f:
            sg = f["summary"]
            for _k in ["nc_delta", "nc_wasserstein", "nc_shift_p", "nc_shift_q",
                       "mde_median", "n_reads"]:
                if _k in sg:
                    _summary[_k] = sg[_k][:]

        def _sum_val(key, v_idx, default):
            if key in _summary and v_idx >= 0:
                return _summary[key][int(v_idx)].item()
            return default

        # Ground truth NC (n_reads may be an attr rather than a dataset)
        _gtnc = f["ground_truth_nc"]
        _gtnc_n = (_gtnc["n_reads"][()] if "n_reads" in _gtnc
                   else _gtnc.attrs.get("n_reads", 0))
        gtnc = {
            "nc_vals": _gtnc["nc_vals"][:].tolist(),
            "nc_fracs": _gtnc["nc_fracs"][:].tolist(),
            "n_reads": int(_gtnc_n),
        }

        # Reference sequence (assumed dataset 'metadata/ref_sequence' or fallback)
        if "ref_sequence" in f:
            ref_seq = f["ref_sequence"][()]
            ref_seq = ref_seq.decode() if isinstance(ref_seq, bytes) else str(ref_seq)
        else:
            ref_seq = "N" * meta["analysis_length"]

        # ── Variants ─────────────────────────────────────────
        variants = []
        if "variants" in f:
            for key in f["variants"]:
                grp = f[f"variants/{key}"]
                v_attrs = _attrs(grp)
                position = int(v_attrs.get("position", 0))
                v_idx = int(v_attrs.get("variant_idx", -1))
                v = {
                    "variant_id": v_attrs["variant_id"],
                    "position": position,
                    "position_1based": position + 1,
                    "ref_base": v_attrs["ref_base"],
                    "alt_base": v_attrs["alt_base"],
                    "change_type": v_attrs.get("change_type", "snv"),
                    "n_reads": int(v_attrs.get("n_reads", v_attrs.get("n_nc_matched", 0))),
                    "best_cluster_p": float(v_attrs.get("best_cluster_p", 1.0)),
                    "variant_fdr_q": float(v_attrs.get("variant_fdr_q", 1.0)),
                    "nc_delta": float(_sum_val("nc_delta", v_idx, v_attrs.get("nc_delta", 0))),
                    "nc_wasserstein": float(_sum_val("nc_wasserstein", v_idx, v_attrs.get("nc_wasserstein", 0))),
                    "nc_shift_p": float(_sum_val("nc_shift_p", v_idx, v_attrs.get("nc_shift_p", 1))),
                    "nc_shift_q": float(_sum_val("nc_shift_q", v_idx, v_attrs.get("nc_shift_q", 1))),
                    "mde_median": float(_sum_val("mde_median", v_idx, v_attrs.get("mde_median", 0))),
                    "per_bin": {},
                }
                best_eff = 0.0
                for b in bins:
                    bg = grp[b]
                    delta_obs = bg["delta_obs"][:]
                    variant_occ = bg["variant_occ"][:]
                    q_values = bg["q_values"][:]
                    # accumulate wt_occ = variant_occ - delta_obs
                    implied_wt = variant_occ - delta_obs
                    if _wt_acc[b] is None:
                        _wt_acc[b] = implied_wt.astype(np.float64)
                    else:
                        _wt_acc[b] += implied_wt
                    pb = {
                        "delta_obs": _round_list(delta_obs.tolist(), precision),
                        "variant_occ": _round_list(variant_occ.tolist(), precision),
                        "empirical_p": _round_list(bg["empirical_p"][:].tolist(), precision),
                        "q_values": _round_list(q_values.tolist(), precision),
                        "mde": _round_list(bg["mde"][:].tolist(), precision),
                        "max_abs_delta": float(np.max(np.abs(delta_obs))),
                        "n_sig_positions": int(np.sum(q_values < 0.10)),
                        "significant_clusters": [],
                    }
                    if "significant_clusters" in bg:
                        for cl_key in bg["significant_clusters"]:
                            cl = bg[f"significant_clusters/{cl_key}"]
                            ca = _attrs(cl)
                            pb["significant_clusters"].append({
                                "bin_label": b,
                                "abs_start": int(ca["abs_start"]),
                                "abs_end": int(ca["abs_end"]),
                                "width": int(ca["width"]),
                                "sum_abs_delta": float(ca.get("sum_abs_delta", 0)),
                                "max_abs_delta": float(ca.get("max_abs_delta", 0)),
                                "mean_signed_delta": float(ca.get("mean_signed_delta", 0)),
                                "direction": ca.get("direction", "loss"),
                                "sign_consistency": float(ca.get("sign_consistency", 0)),
                                "is_motif": bool(ca.get("is_motif", False)),
                                "variant_distance": int(ca.get("variant_distance", 0)),
                                "ref_sequence": ca.get("ref_sequence", ""),
                                "peak_position": int(ca.get("peak_position", ca["abs_start"])),
                            })
                    pb["n_sig_clusters"] = len(pb["significant_clusters"])
                    v["per_bin"][b] = pb
                    if pb["max_abs_delta"] > abs(best_eff):
                        # signed magnitude — pull the sign from the largest |Δ|
                        idx = int(np.argmax(np.abs(delta_obs)))
                        best_eff = float(delta_obs[idx])
                v["best_effect"] = best_eff
                variants.append(v)

        # Finalise reconstructed wt_occupancy (mean across variants)
        if not wt and variants:
            n = len(variants)
            wt = {b: _round_list((_wt_acc[b] / n).tolist(), precision)
                  for b in bins if _wt_acc[b] is not None}

        # ── Motifs ────────────────────────────────────────────
        motifs = []
        density_tracks = {}
        if "motifs" in f:
            mg = f["motifs"]
            if "density_tracks" in mg:
                for key in mg["density_tracks"]:
                    density_tracks[key] = _round_list(
                        mg[f"density_tracks/{key}"][:].tolist(), precision)
            if "calls" in mg:
                for mkey in mg["calls"]:
                    cg = mg[f"calls/{mkey}"]
                    ma = _attrs(cg)
                    sens_count = cg["sensitivity_count"][:].tolist()
                    sens_mean = cg["sensitivity_mean_signed_delta"][:].tolist()
                    contributing = cg["contributing_variant_ids"][:]
                    contributing = [c.decode() if isinstance(c, bytes) else c for c in contributing]
                    motifs.append({
                        "m_id": mkey,
                        "bin": ma["bin"],
                        "direction": ma["direction"],
                        "abs_start": int(ma["abs_start"]),
                        "abs_end": int(ma["abs_end"]),
                        "width": int(ma["width"]),
                        "ref_sequence": ma["ref_sequence"],
                        "n_variants": int(ma["n_variants"]),
                        "peak_density": int(ma["peak_density"]),
                        "base_order": ma.get("base_order", "ACGT"),
                        "sensitivity_count": sens_count,
                        "sensitivity_mean_signed_delta": sens_mean,
                        "contributing_variant_ids": contributing,
                        "tf_label": ma.get("tf_label", mkey),
                        "tf_db": ma.get("tf_db", "—"),
                    })

        # ── Co-occupancy ──────────────────────────────────────
        cooccupancy = []
        if "cooccupancy" in f:
            for pkey in f["cooccupancy"]:
                pg = f[f"cooccupancy/{pkey}"]
                pa = _attrs(pg)
                cooccupancy.append({
                    "pair_id": pkey,
                    "site1_id": pa["site1_id"],
                    "site2_id": pa["site2_id"],
                    "site1_bin": pa["site1_bin"],
                    "site2_bin": pa["site2_bin"],
                    "site1_abs": _to_pylist(pa["site1_abs"]),
                    "site2_abs": _to_pylist(pa["site2_abs"]),
                    "site1_label": pa.get("site1_label", pa["site1_id"]),
                    "site2_label": pa.get("site2_label", pa["site2_id"]),
                    "n_instruments": int(pa["n_instruments"]),
                    "weighted_mean_excess": float(pa["weighted_mean_excess"]),
                    "p_two_sided": float(pa["p_two_sided"]),
                    "fdr_q": float(pa["fdr_q"]),
                    "frac_instruments_consistent": float(pa["frac_instruments_consistent"]),
                    "mde": float(pa["mde"]),
                    "underpowered": bool(pa["underpowered"]),
                    "wt_p2_given_1": float(pa["wt_p2_given_1"]),
                    "wt_p2_given_0": float(pa["wt_p2_given_0"]),
                    "is_call": bool(pa["is_call"]),
                    "validated_call": bool(pa["validated_call"]),
                    "secondary_artifact": bool(pa.get("secondary_artifact", False)),
                    "sec_max_secondary_freq": float(pa.get("sec_max_secondary_freq", 0)),
                })

        # ── TF candidates (TSV companion if present) ──────────
        tf_candidates = []
        tf_tsv = h5_path.with_name(f"{h5_path.stem}_tf_candidates.tsv")
        if tf_tsv.exists():
            import csv
            with open(tf_tsv) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    tf_candidates.append({
                        "motif_idx": int(row["motif_idx"]),
                        "bin": row["bin"],
                        "abs_start": int(row["abs_start"]),
                        "abs_end": int(row["abs_end"]),
                        "motif_seq": row["motif_seq"],
                        "db": row["db"],
                        "tf": row["tf"],
                        "matrix_id": row["matrix_id"],
                        "pwm_len": int(row["pwm_len"]),
                        "strand": row["strand"],
                        "offset": int(row["offset"]),
                        "logodds": float(row["logodds"]),
                        "score_norm": float(row["score_norm"]),
                        "pwm_ic_mean": float(row["pwm_ic_mean"]),
                        "ic_sensitivity_concordance": float(row["ic_sensitivity_concordance"]) if row.get("ic_sensitivity_concordance", "") != "" else None,
                    })

        # ── External MPRA scores (optional) ────────────────────
        mpra_records = []
        if mpra_tsv is not None and Path(mpra_tsv).exists():
            import csv
            with open(mpra_tsv) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    mpra_records.append({
                        "position": int(row["position"]),
                        "ref": row.get("ref", ""),
                        "alt": row.get("alt", ""),
                        "score": float(row.get("score", 0)),
                        "pvalue": float(row.get("pvalue", 1)),
                    })

    return {
        "metadata": meta,
        "ref_sequence": ref_seq,
        "wt_occupancy": wt,
        "ground_truth_nc": gtnc,
        "variants": variants,
        "motifs": motifs,
        "density_tracks": density_tracks,
        "cooccupancy": cooccupancy,
        "tf_candidates": tf_candidates,
        "mpra": mpra_records,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("h5", nargs="+", type=Path, help="One or more {sample}.h5 files")
    ap.add_argument("--mpra", type=Path, help="Optional external MPRA TSV")
    ap.add_argument("-o", "--out", type=Path, default=Path("data.json"))
    ap.add_argument("--precision", type=int, default=4,
                    help="Decimal places for per-position float arrays (default: 4)")
    args = ap.parse_args()

    samples = {}
    order = []
    for p in args.h5:
        s = convert(p, mpra_tsv=args.mpra, precision=args.precision)
        name = s["metadata"]["sample_name"]
        samples[name] = s
        order.append(name)

    out = {"samples": samples, "sampleOrder": order}
    args.out.write_text(json.dumps(out))
    print(f"Wrote {args.out} ({len(samples)} sample(s))", file=sys.stderr)


if __name__ == "__main__":
    main()
