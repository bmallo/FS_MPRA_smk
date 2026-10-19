#!/usr/bin/env python3
"""Bounded, read-only Rep3 caller and m6A provenance audit."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

import pysam


MM_HEADER = re.compile(r"([ACGTUN])([+-])([a-z]+|[0-9]+)([.?]?)\Z")


def ml_groups(read):
    """Consume MM/ML jointly, including non-m6A groups; never use modified_bases."""
    if not read.has_tag("MM") or not read.has_tag("ML"):
        return None
    ml = read.get_tag("ML")
    if getattr(ml, "typecode", None) != "B":
        raise ValueError("ML is not B:C")
    cursor = 0
    target = []
    all_scores = []
    groups = []
    for section in read.get_tag("MM").split(";"):
        if not section:
            continue
        header, *deltas = section.split(",")
        match = MM_HEADER.fullmatch(header)
        if not match:
            raise ValueError(f"Unsupported MM header {header!r}")
        base, strand, codes, omitted = match.groups()
        mods = [codes] if codes.isdigit() else list(codes)
        if any(not delta.isdigit() for delta in deltas):
            raise ValueError("Non-integer MM delta")
        end = cursor + len(deltas) * len(mods)
        if end > len(ml):
            raise ValueError("MM exceeds ML")
        for row in range(len(deltas)):
            for col, code in enumerate(mods):
                score = int(ml[cursor + row * len(mods) + col])
                all_scores.append(score)
                if (base, strand, code) in (("A", "+", "a"), ("T", "-", "a")):
                    target.append(score)
        groups.append({"header": header, "calls": len(deltas), "omitted": omitted})
        cursor = end
    if cursor != len(ml):
        raise ValueError("Extra ML values")
    return target, all_scores, groups


def ma_nuc_aq(read):
    if not read.has_tag("MA"):
        return None
    parts = read.get_tag("MA").split(";")[1:]
    aq = read.get_tag("AQ") if read.has_tag("AQ") else ()
    cursor = 0
    nuc = []
    labels = Counter()
    for part in parts:
        head, sep, data = part.partition(":")
        if not sep:
            raise ValueError("Malformed MA section")
        match = re.fullmatch(r"([^+.-]+)[+.-](.*)", head)
        if not match:
            raise ValueError("Malformed MA header")
        name, spec = match.groups()
        count = len([x for x in data.split(",") if x])
        labels[name] += count
        for _ in range(count):
            if cursor + len(spec) > len(aq):
                raise ValueError("MA requires more AQ bytes")
            if name == "nuc":
                nuc.append(int(aq[cursor]) if spec else None)
            cursor += len(spec)
    if cursor != len(aq):
        raise ValueError("Extra AQ bytes")
    return nuc, labels


def bam_identity(path):
    p = Path(path)
    st = p.stat()
    return {"path": str(p), "device": st.st_dev, "inode": st.st_ino,
            "bytes": st.st_size, "mtime_ns": st.st_mtime_ns, "ctime_ns": st.st_ctime_ns}


def check_bam_identity(path, expected):
    if bam_identity(path) != expected:
        raise RuntimeError(f"BAM changed during audit: {path}")


def scan(path, limit, expected):
    p = Path(path)
    check_bam_identity(p, expected)
    out = {**expected,
           "records": 0, "with_MM_ML": 0, "without_MM_ML": 0, "target_m6a_ML": {},
           "all_ML": {}, "non_m6a_ML": {}, "MM_groups": {}, "with_MA": 0, "MA_calls": {},
           "nuc_AQ": {}, "PG": []}
    target_hist, all_hist, group_hist, ma_hist, nuc_hist = (Counter() for _ in range(5))
    with pysam.AlignmentFile(str(p), "rb", threads=1, check_sq=False) as bam:
        out["PG"] = bam.header.to_dict().get("PG", [])
        for read in bam.fetch(until_eof=True):
            if out["records"] >= limit:
                break
            out["records"] += 1
            result = ml_groups(read)
            if result is None:
                out["without_MM_ML"] += 1
            else:
                out["with_MM_ML"] += 1
                target, all_scores, groups = result
                target_hist.update(target)
                all_hist.update(all_scores)
                group_hist.update((g["header"] for g in groups))
            ma = ma_nuc_aq(read)
            if ma is not None:
                out["with_MA"] += 1
                nuc, labels = ma
                ma_hist.update(labels)
                nuc_hist.update("missing" if x is None else str(x) for x in nuc)
    check_bam_identity(p, expected)
    for key, hist in (("target_m6a_ML", target_hist), ("all_ML", all_hist),
                      ("non_m6a_ML", all_hist - target_hist),
                      ("MM_groups", group_hist), ("MA_calls", ma_hist), ("nuc_AQ", nuc_hist)):
        out[key] = dict(sorted(hist.items(), key=lambda pair: str(pair[0])))
    return out


def small_file(path):
    p = Path(path)
    before = bam_identity(p)
    digest = hashlib.sha256(p.read_bytes()).hexdigest()
    check_bam_identity(p, before)
    return {**before, "sha256": digest}


def check_small_file(path, expected):
    if small_file(path) != expected:
        raise RuntimeError(f"Source file changed during audit: {path}")


def write_atomic_new(path, data):
    """Publish complete output without replacing any existing audit."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Audit output already exists: {path}")
    fd, tmp = tempfile.mkstemp(prefix=".caller_audit_", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(tmp, path)
    finally:
        os.unlink(tmp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=1000)
    ap.add_argument("--output", type=Path, default=Path("results/rep3_caller_provenance_20261018/audit_v2.json"))
    args = ap.parse_args()
    if not 1 <= args.limit <= 1000:
        ap.error("limit must be 1–1000")
    if args.output.exists():
        raise FileExistsError(f"Audit output already exists: {args.output}")
    cfg_path = Path("results/rep3_ledger_20261016/config.json")
    contract_path = Path("config/rep3_existing_footprints.json")
    initial = {str(p): small_file(p) for p in (cfg_path, contract_path)}
    cfg = json.loads(cfg_path.read_text())
    contract = json.loads(contract_path.read_text())
    aligned = Path(cfg["project_root"]) / "fiberseq/HepG2/replicate_3/aligned/PS01898.hg38_LDLR.ft.aligned.bam"
    paths = [cfg["input"]["bam"], str(aligned), contract["bam"]]
    caller_logs = [Path(p) for p in contract["evidence"]["caller_logs"]]
    source_paths = [cfg_path, contract_path,
                    Path("/mmfs1/gscratch/stergachislab/bmallo/large_home/data/pacbio/r84046_20261002_184948_1_A01/logs/ft_r84046_20261002_184948_1_A01.sh"),
                    Path("/mmfs1/home/bmallo/miniconda3/envs/general/lib/python3.12/site-packages/fiberhmm/inference/tf_recaller.py"),
                    Path("/mmfs1/home/bmallo/miniconda3/envs/general/lib/python3.12/site-packages/fiberhmm/inference/fused_stages.py"),
                    Path("/mmfs1/home/bmallo/miniconda3/envs/general/lib/python3.12/site-packages/fiberhmm/cli/call.py"),
                    Path("/mmfs1/home/bmallo/miniconda3/envs/general/lib/python3.12/site-packages/fiberhmm/io/ma_tags.py"),
                    Path(__file__).resolve(), *caller_logs]
    sources = {str(p): initial[str(p)] if str(p) in initial else small_file(p)
               for p in source_paths}
    bams = {str(p): bam_identity(p) for p in paths}
    result = {"schema_version": "fsmpra.rep3_caller_provenance.v2", "limit_per_bam": args.limit,
              "source_check": "BAM device/inode/size/mtime/ctime and small-file SHA256 verified before and after; BAM contents not fully hashed",
              "bam_samples": [scan(path, args.limit, bams[str(path)]) for path in paths],
              "caller_logs": [{**small_file(p),
                               "banner": [line.strip() for line in p.read_text(errors="replace")[:650].splitlines()
                                          if line.strip() and not line.startswith("====")],
                               "mentions_with_scores": "with-scores" in p.read_text(errors="replace")}
                              for p in caller_logs],
              "source_files": [sources[str(p)] for p in source_paths if p not in caller_logs]}
    for path in paths:
        check_bam_identity(path, bams[str(path)])
    for p in source_paths:
        check_small_file(p, sources[str(p)])
    write_atomic_new(args.output, json.dumps(result, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
