"""Auditable raw-BAM processing. Sequence assignment never uses chromatin outcomes.

The processing schema is separate from the older descriptive pilot schema.
All stages publish new directories atomically and preserve original read IDs.
"""

import csv
import hashlib
import json
import os
import re
import sqlite3
import subprocess
import tempfile
from collections import Counter
from functools import lru_cache
from pathlib import Path

from .catalog import build_library
from .features import interval_runs, parse_ma, projected_coverage
from .io import atomic_output, open_text, read_json, source_info, write_json, write_table
from .storage import scratch_database

SCHEMA = "fsmpra.processing.v1"
OBS_FIELDS = ["read_id", "sample_id", "library_id", "observed_sequence", "variant_id",
              "genotype_class", "genotype_status", "barcode", "barcode_qualities",
              "barcode_status", "passes", "read_quality", "mapq", "read_length",
              "strand", "plasmid_class", "circular", "m6a_status", "at_count", "m6a_count",
              "methylation_fraction", "sequence_qc", "length_qc", "methylation_qc",
              "analysis_eligible"]


def config_read(path):
    cfg = read_json(path)
    if cfg.get("schema_version") != SCHEMA:
        raise ValueError("Unsupported processing configuration")
    lib = cfg["library"]
    a, b = lib["barcode_start0"], lib["barcode_end0"]
    if not 0 < a < b or b - a != lib["barcode_length"]:
        raise ValueError("Invalid barcode coordinates/length")
    if lib["barcode_orientation"] not in ("reference", "reverse_complement"):
        raise ValueError("Set an explicit barcode orientation")
    q = cfg["qc"]
    if (q["min_passes"] < 0 or not 0 <= q["min_read_quality"] <= 1 or
            not 0 <= q["min_mapq"] <= 255 or not 0 <= q["min_base_quality"] <= 93 or
            not 0 <= q["ml_threshold"] <= 255 or q["end_slop"] < 0 or
            q["max_softclip"] < 0 or q["concatemer_length_fraction"] <= 1):
        raise ValueError("Invalid QC thresholds")
    if not 0 <= q["methylation_min"] <= q["methylation_max"] <= 1:
        raise ValueError("Invalid methylation interval")
    if cfg["input"]["modifications"] not in ("existing", "predict"):
        raise ValueError("modifications must be existing or predict")
    return cfg


def check_sources(sources):
    for original in sources:
        if source_info(original["path"], "sha256" in original) != original:
            raise ValueError("Input changed during processing: " + original["path"])


def env_for(scratch):
    env = dict(os.environ)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
        env[key] = "1"
    env["NUMBA_CACHE_DIR"] = str(Path(scratch) / "numba")
    env["TMPDIR"] = str(scratch)
    return env


def run_command(argv, log, scratch):
    """No shell interpolation; tools get explicit resources and disposable cache."""
    with open(log, "a") as handle:
        handle.write(json.dumps([str(a) for a in argv]) + "\n")
        handle.flush()
        result = subprocess.run([str(a) for a in argv], stdout=handle, stderr=handle,
                                env=env_for(scratch))
    if result.returncode:
        with open(log, "rb") as handle:
            handle.seek(max(0, os.path.getsize(log) - 6000))
            tail = handle.read().decode(errors="replace")
        raise RuntimeError(f"External command exited {result.returncode}:\n{tail}")


def tool_identity(cfg, name):
    tool = cfg["tools"][name]
    if name == "pbmm2":
        version = subprocess.check_output([tool, "--version"], text=True).splitlines()[0].split()[-1]
    elif name == "fiberhmm":
        # Query metadata through the installed entrypoint's interpreter, without
        # importing the caller (which would initialize numba caches).
        with open(tool) as handle:
            interpreter = handle.readline().strip().removeprefix("#!")
        if not Path(interpreter).is_file():
            raise ValueError("FiberHMM entrypoint needs an absolute Python shebang")
        version = subprocess.check_output([interpreter, "-c", "import importlib.metadata; print(importlib.metadata.version('fiberhmm'))"], text=True).strip()
    else:
        raise ValueError("Unsupported version check")
    if version != cfg["tool_versions"][name]:
        raise ValueError(f"{name} version {version} differs from pinned {cfg['tool_versions'][name]}")
    return {"version": version, "entrypoint": source_info(tool, True)}


def preflight(config_path, output, inspect_reads=1000):
    import pysam
    cfg = config_read(config_path)
    if inspect_reads < 1:
        raise ValueError("inspect_reads must be positive")
    sources = [source_info(config_path, True), source_info(cfg["input"]["bam"])]
    seen, calls, kinetics = 0, 0, 0
    with pysam.AlignmentFile(cfg["input"]["bam"], "rb", check_sq=False) as bam:
        programs = bam.header.to_dict().get("PG", [])
        for read in bam:
            if seen == inspect_reads:
                break
            if not read.is_unmapped or read.is_secondary or read.is_supplementary:
                raise ValueError("This entrypoint requires unaligned primary HiFi reads")
            seen += 1
            mm = read.get_tag("MM") if read.has_tag("MM") else ""
            calls += int(read.has_tag("ML") and any(s.startswith(("A+a", "T-a")) for s in mm.split(";")))
            kinetics += int(all(read.has_tag(t) for t in ("fi", "fp", "ri", "rp")))
    if not seen or (cfg["input"]["modifications"] == "existing" and not calls):
        raise ValueError("Missing m6A MM/ML calls in inspected input; select prediction only if kinetics exist")
    if cfg["input"]["modifications"] == "predict" and kinetics != seen:
        raise ValueError("m6A prediction requires fi/fp/ri/rp kinetics")
    with atomic_output(output) as out:
        report = {"schema_version": SCHEMA, "sources": sources, "inspected_reads": seen,
                  "reads_with_calls": calls, "reads_with_kinetics": kinetics, "programs": programs,
                  "scope": "bounded preflight; every plasmid record is checked again during extraction"}
        check_sources(sources)
        write_json(out / "report.json", report)
    return report


def prepare(config_path, output):
    import pysam
    cfg = config_read(config_path)
    lib = build_library({"project_root": cfg["project_root"]}, cfg["library"])
    # FastaFile uses the existing index; it never loads the human genome into RAM.
    fasta = cfg["alignment"]["reference_fasta"]
    with pysam.FastaFile(fasta) as ref:
        seq = ref.fetch(lib["contig"]).upper()
    if hashlib.sha256(seq.encode()).hexdigest() != lib["reference_sequence_sha256"]:
        raise ValueError("Combined alignment reference disagrees with library")
    if lib["barcode_end0"] >= len(seq):
        raise ValueError("Barcode needs flanking sequence on both sides")
    if set(seq[lib["barcode_start0"]:lib["barcode_end0"]]) != {"N"}:
        raise ValueError("Expected N-masked barcode interval in plasmid reference")
    if not (lib["mutable_ref_end0"] <= lib["barcode_start0"] or
            lib["barcode_end0"] <= lib["mutable_ref_start0"]):
        raise ValueError("Mutable element and barcode overlap")
    sources = [source_info(config_path, True), *lib["sources"], source_info(fasta), source_info(fasta + ".fai")]
    with atomic_output(output) as out:
        write_json(out / "catalog.json", {"schema_version": SCHEMA, "library": lib, "sources": sources})
        with open(out / "design.fa", "w") as handle:
            for v in lib["variants"]:
                handle.write(f">{v['variant_id']}\n{v['sequence']}\n")
        write_table(out / "oligo_aliases.tsv", ({"element_id": v["variant_id"], "oligo_name": name}
                    for v in lib["variants"] for name in v["oligo_names"]), ["element_id", "oligo_name"])
        write_table(out / "elements.tsv", ({k: v[k] for k in ("variant_id", "class", "source_variant", "sequence", "sequence_sha256")}
                    for v in lib["variants"]), ["variant_id", "class", "source_variant", "sequence", "sequence_sha256"])
        (out / "plasmid.fa").write_text(f">{lib['contig']}\n{seq}\n")
        check_sources(sources)
    return {"elements": len(lib["variants"]), "output": str(output)}


def align(config_path, output, threads=1):
    import pysam
    cfg = config_read(config_path)
    if threads < 1:
        raise ValueError("threads must be positive")
    identity = tool_identity(cfg, "pbmm2")
    alignment = cfg["alignment"]
    reference = alignment.get("index") or alignment["reference_fasta"]
    sources = [source_info(config_path, True), source_info(cfg["input"]["bam"]),
               source_info(reference), source_info(alignment["reference_fasta"])]
    with atomic_output(output) as out, tempfile.TemporaryDirectory(dir=cfg["scratch_dir"], prefix="fsmpra-align-") as tmp:
        input_bam = cfg["input"]["bam"]
        if cfg["input"]["modifications"] == "predict":
            predicted = str(Path(tmp) / "modified.bam")
            run_command([cfg["tools"]["ft"], "predict-m6a", "-t", threads, input_bam, predicted], out / "commands.log", tmp)
            input_bam = predicted
        unsorted = str(Path(tmp) / "aligned.bam")
        run_command([cfg["tools"]["pbmm2"], "align", "--preset", "HIFI", "--unmapped", "-j", threads,
                     reference, input_bam, unsorted], out / "commands.log", tmp)
        # Sequential sort avoids adding sort threads on top of alignment workers.
        pysam.sort("-@", "0", "-m", "1G", "-T", str(Path(tmp) / "sort"), "-o", str(out / "aligned.bam"), unsorted)
        pysam.index(str(out / "aligned.bam"))
        with pysam.AlignmentFile(out / "aligned.bam", "rb") as bam, pysam.FastaFile(alignment["reference_fasta"]) as ref:
            if dict(zip(bam.references, bam.lengths)) != dict(zip(ref.references, ref.lengths)):
                raise ValueError("Alignment index and configured reference dictionaries disagree")
            programs = bam.header.to_dict().get("PG", [])
        check_sources(sources)
        write_json(out / "audit.json", {"schema_version": SCHEMA, "sources": sources, "threads": threads,
                   "programs": programs, "tool": identity,
                   "index_identity": "path/size/mtime plus reference dictionary; index sequence identity is not proven"})
    return {"output": str(output)}


def extract_interval(read, start, end, reference, flank=4):
    """Slice between aligned outside anchors, retaining insertions and deletions.

    BAM stored SEQ is already in reference orientation, even for reverse reads.
    Require exact short flanks; uncertain boundaries are not silently repaired.
    """
    if read.query_sequence is None or start < flank or end + flank > len(reference):
        return "", [], "missing_sequence_or_flank"
    # Walk CIGAR blocks instead of allocating a tuple for every base in every
    # plasmid. Only the eight flanking aligned positions are needed here.
    anchors = {}
    query_pos, ref_pos = 0, read.reference_start
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):  # M, =, X consume query and reference
            for lo, hi in ((start - flank, start), (end, end + flank)):
                for r in range(max(ref_pos, lo), min(ref_pos + length, hi)):
                    anchors[r] = query_pos + r - ref_pos
            query_pos += length
            ref_pos += length
        elif op in (1, 4):  # insertion, soft clip
            query_pos += length
        elif op in (2, 3):  # deletion, reference skip
            ref_pos += length
        elif op not in (5, 6):  # hard clip and padding consume neither
            raise ValueError("Unsupported CIGAR operation")
        if ref_pos >= end + flank:
            break
    left = [anchors.get(r) for r in range(start - flank, start)]
    right = [anchors.get(r) for r in range(end, end + flank)]
    if any(q is None for q in left + right):
        return "", [], "incomplete_flanks"
    seq = read.query_sequence
    if (left != list(range(left[0], left[0] + flank)) or
            right != list(range(right[0], right[0] + flank)) or
            seq[left[0]:left[-1] + 1] != reference[start - flank:start] or
            seq[right[0]:right[-1] + 1] != reference[end:end + flank]):
        return "", [], "flank_disagreement"
    a, b = left[-1] + 1, right[0]
    qualities = list(read.query_qualities[a:b]) if read.query_qualities is not None else []
    return seq[a:b], qualities, "ok"


@lru_cache(maxsize=256)
def _probability_table(threshold):
    return bytes(int(i >= threshold) for i in range(256))


def methylation(read, threshold):
    """Count m6A without expanding every modification into a genomic position.

    Nonnegative MM skips are strictly increasing within a group, so checking the
    last opportunity validates every index. A and T groups cannot overlap.
    Repeated target groups use the reference decoder to deduplicate positions.
    """
    seq = read.query_sequence or ""
    at = seq.count("A") + seq.count("T")
    if not read.has_tag("MM") or not read.has_tag("ML"):
        return "missing", at, None, None
    mm = read.get_tag("MM")
    sections = [s.split(",", 1)[0] for s in mm.split(";") if s.startswith(("A+a", "T-a"))]
    if not sections:
        return "missing_m6a", at, None, None
    if read.has_tag("MN") and read.get_tag("MN") != len(seq):
        raise ValueError("MN/SEQ length mismatch: " + read.query_name)
    ml_array = read.get_tag("ML")
    if ml_array.typecode != "B":
        return _methylation_reference(read, threshold)
    ml = ml_array.tobytes()
    forward = read.get_forward_sequence()
    cursor, calls = 0, 0
    targets = set()
    for section in mm.split(";"):
        if not section:
            continue
        header, *encoded = section.split(",")
        match = re.fullmatch(r"([ACGTUN])([+-])([a-z]+|[0-9]+)([.?]?)", header)
        if not match:
            raise ValueError("Unsupported MM group: " + header)
        base, strand, codes, _ = match.groups()
        mods = [codes] if codes.isdigit() else list(codes)
        deltas = list(map(int, encoded))
        n = len(deltas)
        if n and min(deltas) < 0:
            raise ValueError("Negative MM delta")
        opportunities = len(forward) if base == "N" else forward.count(base)
        end = cursor + n * len(mods)
        if sum(deltas) + n > opportunities or end > len(ml):
            raise ValueError("MM/ML/SEQ length disagreement: " + read.query_name)
        for j, code in enumerate(mods):
            if (base, strand, code) in (("A", "+", "a"), ("T", "-", "a")):
                if base in targets:
                    return _methylation_reference(read, threshold)
                targets.add(base)
                calls += ml[cursor+j:end:len(mods)].translate(_probability_table(threshold)).count(b"\x01")
        cursor = end
    if cursor != len(ml):
        raise ValueError("Extra ML probabilities: " + read.query_name)
    if any(s.endswith("?") for s in sections) or not at:
        return "unknown_opportunities", at, calls, None
    return "ok", at, calls, calls / at


def _methylation_reference(read, threshold):
    seq = read.query_sequence or ""
    at = seq.count("A") + seq.count("T")
    if not read.has_tag("MM") or not read.has_tag("ML"):
        return "missing", at, None, None
    sections = [s.split(",")[0] for s in read.get_tag("MM").split(";") if s.startswith(("A+a", "T-a"))]
    if not sections:
        return "missing_m6a", at, None, None
    # Decode the two m6A channels explicitly: the system pysam/htslib's
    # modified_bases accessor segfaults on real mixed m6A/5mC/5hmC records.
    # MM deltas use original forward SEQ, regardless of alignment orientation.
    if read.has_tag("MN") and read.get_tag("MN") != len(seq):
        raise ValueError("MN/SEQ length mismatch: " + read.query_name)
    forward = read.get_forward_sequence()
    ml = read.get_tag("ML")
    cursor = 0
    positions = set()
    for section in read.get_tag("MM").split(";"):
        if not section:
            continue
        header, *deltas = section.split(",")
        match = re.fullmatch(r"([ACGTUN])([+-])([a-z]+|[0-9]+)([.?]?)", header)
        if not match:
            raise ValueError("Unsupported MM group: " + header)
        base, strand, codes, _ = match.groups()
        mods = [codes] if codes.isdigit() else list(codes)
        opportunities = [i for i, letter in enumerate(forward) if letter == base or base == "N"]
        index = -1
        for delta in deltas:
            skip = int(delta)
            if skip < 0:
                raise ValueError("Negative MM delta")
            index += skip + 1
            if index >= len(opportunities) or cursor + len(mods) > len(ml):
                raise ValueError("MM/ML/SEQ length disagreement: " + read.query_name)
            for j, code in enumerate(mods):
                if (base, strand, code) in (("A", "+", "a"), ("T", "-", "a")) and ml[cursor + j] >= threshold:
                    positions.add(opportunities[index])
            cursor += len(mods)
    if cursor != len(ml):
        raise ValueError("Extra ML probabilities: " + read.query_name)
    # '?' denotes unknown omitted positions: never turn those into unmethylated A/T.
    if any(s.endswith("?") for s in sections) or not at:
        return "unknown_opportunities", at, len(positions), None
    return "ok", at, len(positions), len(positions) / at


def structure(read, reference_length, qc):
    softclip = sum(n for op, n in read.cigartuples or [] if op == 4)
    full = read.reference_start <= qc["end_slop"] and read.reference_end >= reference_length - qc["end_slop"]
    pc = "concatemer" if read.query_length > reference_length * qc["concatemer_length_fraction"] else "full" if full else "partial"
    circular = pc == "full" and softclip <= qc["max_softclip"] and not read.has_tag("SA")
    return pc, circular


def partition(config_path, catalog_path, bam_path, output):
    """Caller partition is independent of sequence/barcode/analysis filtering."""
    import pysam
    cfg = config_read(config_path)
    lib = read_json(catalog_path)["library"]
    sources = [source_info(config_path, True), source_info(catalog_path, True), source_info(bam_path)]
    counts = Counter()
    with atomic_output(output) as out, pysam.AlignmentFile(bam_path, "rb") as bam, \
            pysam.AlignmentFile(out / "circular.bam", "wb", template=bam) as circ, \
            pysam.AlignmentFile(out / "linear.bam", "wb", template=bam) as linear, \
            pysam.AlignmentFile(out / "uncallable.bam", "wb", template=bam) as uncallable:
        for read in bam.fetch(lib["contig"]):
            if read.is_secondary or read.is_supplementary:
                continue
            counts["plasmid_primary"] += 1
            _, circular = structure(read, lib["reference_length"], cfg["qc"])
            ms, _, _, _ = methylation(read, cfg["qc"]["ml_threshold"])
            for tag in ("MA", "AQ", "AN", "ns", "nl", "as", "al", "nq", "aq"):
                if read.has_tag(tag):
                    read.set_tag(tag, None)
            mode = "uncallable" if ms != "ok" else "circular" if circular else "linear"
            counts[mode] += 1
            {"uncallable": uncallable, "circular": circ, "linear": linear}[mode].write(read)
        check_sources(sources)
        write_json(out / "qc.json", {"schema_version": SCHEMA, "sources": sources, "counts": dict(counts)})
    return {"output": str(output), "counts": dict(counts)}


def extract(config_path, catalog_path, bam_path, output, write_header=True):
    import pysam
    cfg = config_read(config_path)
    lib = read_json(catalog_path)["library"]
    from .io import read_fasta
    reference = read_fasta(Path(catalog_path).parent / "plasmid.fa", lib["contig"])
    variants = {v["sequence"]: v for v in lib["variants"]}
    q = cfg["qc"]
    sources = [source_info(config_path, True), source_info(catalog_path, True), source_info(bam_path)]
    counts = Counter()
    with atomic_output(output) as out, scratch_database(cfg["scratch_dir"]) as db, \
            pysam.AlignmentFile(bam_path, "rb") as bam:
        db.execute("CREATE TABLE ids (id TEXT PRIMARY KEY) WITHOUT ROWID")
        if bam.get_reference_length(lib["contig"]) != len(reference):
            raise ValueError("BAM reference length mismatch")
        with open_text(out / "observations.tsv.gz", "wt") as handle:
            writer = csv.DictWriter(handle, fieldnames=OBS_FIELDS, delimiter="\t")
            if write_header:
                writer.writeheader()
            for read in bam:
                counts["alignment_records"] += 1
                if read.is_secondary or read.is_supplementary:
                    counts["secondary_or_supplementary"] += 1
                    continue
                # All primary IDs, including genomic/unmapped, must be unique.
                try:
                    db.execute("INSERT INTO ids VALUES (?)", (read.query_name,))
                except sqlite3.IntegrityError as exc:
                    raise ValueError("Duplicate primary read ID: " + read.query_name) from exc
                if counts["alignment_records"] % 10000 == 0:
                    db.commit()
                if read.is_unmapped or read.reference_name != lib["contig"]:
                    counts["unmapped" if read.is_unmapped else "other_contig"] += 1
                    continue
                counts["plasmid_primary"] += 1
                sequence, qualities, gs = extract_interval(read, lib["mutable_ref_start0"], lib["mutable_ref_end0"], reference)
                v = variants.get(sequence)
                genotype_ok = gs == "ok" and v is not None and bool(qualities) and min(qualities) >= q["min_base_quality"]
                genotype_status = (gs if gs != "ok" else "off_design" if v is None else
                                   "low_base_quality" if not genotype_ok else "exact_design")
                barcode, bq, bs = extract_interval(read, lib["barcode_start0"], lib["barcode_end0"], reference)
                if lib["barcode_orientation"] == "reverse_complement":
                    barcode, bq = barcode.translate(str.maketrans("ACGTN", "TGCAN"))[::-1], bq[::-1]
                if bs == "ok":
                    bs = ("invalid_bases" if set(barcode) - set("ACGT") else
                          "length_mismatch" if len(barcode) != lib["barcode_length"] else
                          "low_base_quality" if not bq or min(bq) < q["min_base_quality"] else "ok")
                passes = read.get_tag("np") if read.has_tag("np") else None
                rq = read.get_tag("rq") if read.has_tag("rq") else None
                sequence_ok = passes is not None and passes >= q["min_passes"] and rq is not None and rq >= q["min_read_quality"] and read.mapping_quality >= q["min_mapq"]
                pc, circular = structure(read, len(reference), q)
                ms, at, mc, mf = methylation(read, q["ml_threshold"])
                methylation_ok = ms == "ok" and q["methylation_min"] <= mf <= q["methylation_max"]
                length_ok = circular
                eligible = genotype_ok and sequence_ok and length_ok and methylation_ok and v["class"] in ("WT", "SNV")
                row = dict(zip(OBS_FIELDS, [read.query_name, cfg["sample"]["id"], lib["id"], sequence,
                    v["variant_id"] if v else "", v["class"] if v else "unassigned", genotype_status,
                    barcode, ",".join(map(str, bq)), bs, passes, rq, read.mapping_quality, read.query_length,
                    "-" if read.is_reverse else "+", pc, int(circular), ms, at, mc, mf,
                    int(sequence_ok), int(length_ok), int(methylation_ok), int(eligible)]))
                writer.writerow(row)
                counts["genotype_" + genotype_status] += 1
                counts["barcode_" + bs] += 1
                counts["structure_" + pc] += 1
                counts["analysis_eligible"] += int(eligible)
                counts["uncallable_modifications"] += int(ms != "ok")
        check_sources(sources)
        report = {"schema_version": SCHEMA, "sources": sources, "counts": dict(counts), "qc": q,
                  "table_header": write_header,
                  "barcode_orientation": lib["barcode_orientation"],
                  "genotype_method": "exact observed mutable sequence; no barcode consensus rescue",
                  "methylation_role": "eligibility flag only; all primary plasmids retained",
                  "circular_role": "whole-plasmid boundary model, independent of genotype and methylation"}
        write_json(out / "qc.json", report)
    return report


def footprints(config_path, extracted_dir, output, threads=1):
    import pysam
    cfg = config_read(config_path)
    if threads < 2:
        raise ValueError("FiberHMM requires at least 2 allocated CPUs (worker plus driver); use --cores >=2")
    fp = cfg["footprints"]
    identity = tool_identity(cfg, "fiberhmm")
    model = source_info(fp["model"], True)
    if model["sha256"] != fp["model_sha256"]:
        raise ValueError("Footprint model checksum mismatch")
    sources = [source_info(config_path, True), model] + [source_info(Path(extracted_dir) / (m + ".bam")) for m in ("circular", "linear", "uncallable")]
    counts = Counter()
    with atomic_output(output) as out, tempfile.TemporaryDirectory(dir=cfg["scratch_dir"], prefix="fsmpra-footprints-") as tmp:
        files = []
        for mode in ("circular", "linear"):
            src = Path(extracted_dir) / (mode + ".bam")
            dest = Path(tmp) / (mode + ".bam")
            with pysam.AlignmentFile(src, "rb") as bam:
                first = next(bam, None)
                if first is None:
                    with pysam.AlignmentFile(dest, "wb", template=bam):
                        pass
                else:
                    run_command([cfg["tools"]["fiberhmm"], "-i", src, "-o", dest,
                                 "--enzyme", "hia5", "--seq", "pacbio", "--model", fp["model"],
                                 "--recall-model", fp["model"], "--primary", "--min-mapq", "0",
                                 "--prob-threshold", cfg["qc"]["ml_threshold"], "--min-read-length", "1",
                                 "--min-llr", fp["min_llr"], "--min-opps", fp["min_opps"],
                                 "--unify-threshold", "90", "--emission-uplift", "1",
                                 "--io-threads", "0", "-c", max(1, threads - 2),
                                 *(["--circular"] if mode == "circular" else [])], out / "commands.log", tmp)
            files.append(str(dest))
        files.append(str(Path(extracted_dir) / "uncallable.bam"))
        # Sequential concatenation and sorting; no concurrent I/O thread pools.
        joined = str(Path(tmp) / "joined.bam")
        with pysam.AlignmentFile(files[0], "rb") as template, pysam.AlignmentFile(joined, "wb", template=template) as dest:
            for name in files:
                with pysam.AlignmentFile(name, "rb") as bam:
                    for read in bam:
                        counts["output_reads"] += 1
                        counts["with_MA"] += int(read.has_tag("MA"))
                        dest.write(read)
        pysam.sort("-@", "0", "-m", "1G", "-T", str(Path(tmp) / "sort"), "-o", str(out / "footprints.bam"), joined)
        pysam.index(str(out / "footprints.bam"))
        expected = read_json(Path(extracted_dir) / "qc.json")["counts"].get("plasmid_primary", 0)
        if counts["output_reads"] != expected:
            raise ValueError("Footprint caller changed molecule count")
        check_sources(sources)
        write_json(out / "audit.json", {"schema_version": SCHEMA, "sources": sources, "counts": dict(counts),
                   "caller": identity,
                   "evaluated_tracks": ["tf", "nuc", "msp"], "format": "fiberhmm_ma_sparse"})
    return {"output": str(output), "counts": dict(counts)}


def feature_table(config_path, catalog_path, extracted_dir, barcode_dir, footprint_dir, output):
    """One record per plasmid, including missing calls; coordinates are explicit."""
    import pysam
    cfg = config_read(config_path)
    lib = read_json(catalog_path)["library"]
    start, end = lib["mutable_ref_start0"], lib["mutable_ref_end0"]
    paths = [Path(extracted_dir) / "observations.tsv.gz", Path(barcode_dir) / "molecule_assignments.tsv.gz",
             Path(footprint_dir) / "footprints.bam", Path(catalog_path)]
    sources = [source_info(p) for p in paths]
    with atomic_output(output) as out, scratch_database(cfg["scratch_dir"]) as db:
        db.execute("CREATE TABLE molecules (id TEXT PRIMARY KEY, payload TEXT, barcode TEXT, status TEXT, seen INT DEFAULT 0)")
        with open_text(paths[0]) as handle:
            for i, row in enumerate(csv.DictReader(handle, delimiter="\t"), 1):
                db.execute("INSERT INTO molecules(id,payload) VALUES (?,?)", (row["read_id"], json.dumps(row)))
                if i % 10000 == 0:
                    db.commit()
        with open_text(paths[1]) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                cur = db.execute("UPDATE molecules SET barcode=?,status=? WHERE id=? AND status IS NULL",
                                 (row["corrected_barcode"], row["assignment_status"], row["read_id"]))
                if cur.rowcount != 1:
                    raise ValueError("Invalid barcode assignment read ID")
        db.commit()
        fields = OBS_FIELDS + ["corrected_barcode", "assignment_status", "footprint_status", "eligible_ref_positions0",
                               "tf_small_ref_intervals0", "nucleosome_ref_intervals0", "msp_ref_intervals0", "msp_fraction"]
        counts = Counter()
        with open_text(out / "molecule_features.tsv.gz", "wt") as handle, pysam.AlignmentFile(paths[2], "rb") as bam:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for read in bam:
                found = db.execute("SELECT payload,barcode,status,seen FROM molecules WHERE id=?", (read.query_name,)).fetchone()
                if found is None or found[3] or found[2] is None:
                    raise ValueError("Footprint/observation/assignment identity mismatch")
                db.execute("UPDATE molecules SET seen=1 WHERE id=?", (read.query_name,))
                row = json.loads(found[0])
                row.update(corrected_barcode=found[1], assignment_status=found[2], footprint_status="missing")
                qmap = read.get_reference_positions(full_length=True)
                aligned = {p for p in qmap if p is not None and start <= p < end}
                row["eligible_ref_positions0"] = json.dumps(sorted(aligned), separators=(",", ":"))
                if read.has_tag("MA"):
                    tracks = parse_ma(read.get_tag("MA"), read.query_length)
                    row["footprint_status"] = "evaluated"
                    for name, track, lo, hi in (("tf_small", "tf", 10, 60), ("nucleosome", "nuc", 85, None), ("msp", "msp", 1, None)):
                        covered, _ = projected_coverage(tracks.get(track, []), qmap, start, end, lo, hi)
                        row[name + "_ref_intervals0"] = json.dumps([[a + start, b + start] for a, b in interval_runs(covered)], separators=(",", ":"))
                        if name == "msp" and len(aligned) == end - start:
                            row["msp_fraction"] = len(covered) / (end - start)
                counts[row["footprint_status"]] += 1
                writer.writerow(row)
                if sum(counts.values()) % 10000 == 0:
                    db.commit()
        if db.execute("SELECT count(*) FROM molecules WHERE seen=0").fetchone()[0]:
            raise ValueError("Footprint output omitted molecules")
        check_sources(sources)
        write_json(out / "audit.json", {"schema_version": SCHEMA, "sources": sources, "counts": dict(counts),
                   "coordinates": "reference 0-based half-open", "interval": [start, end],
                   "missing_tracks": "Absent section in present MA is evaluated zero; absent MA is missing",
                   "uncertainty": "descriptive features only; no statistical inference"})
    return {"output": str(output), "counts": dict(counts)}
