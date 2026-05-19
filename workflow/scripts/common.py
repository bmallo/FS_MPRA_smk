#!/usr/bin/env python3
"""
common.py — Shared utilities for the FS-MPRA pipeline.

Contains constants, BAM tag parsing, the ReadData class, BAM parsing,
variant grouping, NC-matched subsampling, shared memory helpers,
cluster detection, null calibration, and coverage grid utilities.

These functions are shared between 03_analyze_library.py and
extras/mpra_analysis_qc.py.
"""

import array as arr_module
import hashlib
import json
import logging
import re
import struct
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import shared_memory

import numpy as np
import pysam
from scipy import ndimage
from scipy.stats import norm, wasserstein_distance

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(iterable, **kwargs):
        return iterable

VERSION = "0.2.0"

# ============================================================================
# Constants
# ============================================================================

FOOTPRINT_BINS = {
    'sub_TF': (10, 19, 'sub-TF_10-19bp'),
    'TF':     (20, 40, 'TF_20-40bp'),
    'PIC':    (41, 80, 'PIC_41-80bp'),
    'NUC':    (81, np.inf, 'Nucleosome_81plusbp'),
}

# Track-aware analysis bins (FiberHMM). Each entry:
#   (label, source_track, min_len, max_len)   max_len inclusive.
# tf-track footprints are length sub-binned; the nuc track is a single
# bin keyed by a minimum-length threshold. All bounds are runtime knobs
# (see 03_analyze_library.py); these are the defaults.
NUC_MIN_LEN_DEFAULT = 80
DEFAULT_ANALYSIS_BINS = [
    ('sub_TF', 'tf',  10, 19),
    ('TF',     'tf',  20, 39),
    ('PIC',    'tf',  40, 60),
    ('nuc',    'nuc', NUC_MIN_LEN_DEFAULT, np.inf),
]


def build_analysis_bins(tf_subTF=(10, 19), tf_TF=(20, 39),
                        tf_PIC=(40, 60), nuc_min_len=NUC_MIN_LEN_DEFAULT):
    """Construct the (label, track, min, max) analysis-bin spec from
    configurable bounds. Returns a list parallel to bin_labels."""
    return [
        ('sub_TF', 'tf',  tf_subTF[0], tf_subTF[1]),
        ('TF',     'tf',  tf_TF[0],    tf_TF[1]),
        ('PIC',    'tf',  tf_PIC[0],   tf_PIC[1]),
        ('nuc',    'nuc', nuc_min_len, np.inf),
    ]

CLUSTER_WIDTH_RANGES = {
    'sub-TF_10-19bp':      (5, 30),
    'TF_20-40bp':          (10, 60),
    'PIC_41-80bp':         (25, 100),
    'Nucleosome_81plusbp':  (50, 200),
}

TAG_NUC_COUNT = 'nc'
TAG_PROMOTER_VARIANT = 'PV'
TAG_VARIANT_COUNT = 'VC'
TAG_RAW_VARIANT = 'PR'   # per-read raw (pre-cluster-consensus) calls
TAG_FP_STARTS = 'ns'
TAG_FP_LENGTHS = 'nl'
TAG_FP_QUAL = 'nq'

# FiberHMM structured footprint tag (new format). MA:Z =
#   "<qlen>;nuc+Q:s-l,...;msp+:s-l,...;tf+QQQ:s-l,..."
# starts are 1-based query coords; we expose them 0-based to match the
# legacy `ns` convention used by build_q2r_array/convert_footprints*.
TAG_MA = 'MA'
TAG_MA_QUAL = 'AQ'
# Legacy fibertools binary mirrors, mapped to the correct tracks:
#   ns/nl = nucleosomes, as/al = MSPs (there is no legacy `tf` track).
LEGACY_TRACK_TAGS = {'nuc': ('ns', 'nl'), 'msp': ('as', 'al')}
FIBERHMM_TRACKS = ('nuc', 'msp', 'tf')


# ============================================================================
# Utilities
# ============================================================================

def setup_logging(verbose=False, quiet=False):
    if quiet:
        level = logging.WARNING
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')


def parse_region(region_str):
    """Parse a region string into (chrom, start, end).

    Accepts:
        'CHROM:START-END'  -> (CHROM, START, END)  (1-based inclusive)
        'START-END'        -> (None, START, END)    (1-based inclusive)
        None               -> (None, None, None)

    Returns (chrom_or_None, start_int_or_None, end_int_or_None).
    """
    if region_str is None:
        return None, None, None
    region_str = region_str.strip()
    if ':' in region_str:
        chrom, coords = region_str.split(':', 1)
        parts = coords.split('-')
        if len(parts) != 2:
            raise ValueError(f"Invalid region format: {region_str}. "
                             f"Expected CHROM:START-END or START-END.")
        return chrom, int(parts[0]), int(parts[1])
    elif '-' in region_str:
        parts = region_str.split('-')
        if len(parts) != 2:
            raise ValueError(f"Invalid region format: {region_str}. "
                             f"Expected CHROM:START-END or START-END.")
        return None, int(parts[0]), int(parts[1])
    else:
        # Bare chromosome name
        return region_str, None, None


def safe_hdf5_name(name):
    return (name.replace(' ', '_').replace('(', '').replace(')', '')
                .replace('+', 'plus').replace('>', 'to').replace(':', '_')
                .replace('[', '').replace(']', '').replace('"', '')
                .replace(',', '_').replace('/', '_'))


def safe_variant_hdf5_name(vid):
    """Make a variant ID safe for HDF5 group names."""
    return (vid.replace(':', '_').replace('>', 'to').replace('+', 'plus')
               .replace(' ', '').replace('"', '').replace('[', '')
               .replace(']', '').replace(',', '_').replace('/', '_'))


def benjamini_hochberg(p_values):
    p = np.asarray(p_values, dtype=np.float64)
    n = len(p)
    if n == 0:
        return np.array([], dtype=np.float64)
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n + 1)
    q = p * n / ranks
    sorted_q = q[order]
    for i in range(n - 2, -1, -1):
        sorted_q[i] = min(sorted_q[i], sorted_q[i + 1])
    q[order] = sorted_q
    return np.minimum(q, 1.0)


def derive_sample_name(bam_path):
    """Derive a clean sample name from a BAM file path."""
    import os
    name = os.path.splitext(os.path.basename(bam_path))[0]
    for suffix in ['.sorted', '.aligned', '.merged', '.tagged',
                   '_sorted', '_aligned', '_merged', '_tagged',
                   '.plasmid', '_plasmid']:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def stable_variant_seed(base_seed, vid):
    """Deterministic per-variant seed derived from base_seed and the variant id.

    Python's builtin hash() is salted per process (PYTHONHASHSEED), so
    `base_seed + hash(vid)` is not reproducible across runs and makes the
    serial and parallel paths diverge. Use a stable digest instead.
    """
    digest = hashlib.blake2b(str(vid).encode(), digest_size=8).digest()
    return (int(base_seed) + int.from_bytes(digest, 'little')) % (2**31)


# ============================================================================
# BAM Tag Parsing
# ============================================================================

def parse_tag_array(tag_value):
    if isinstance(tag_value, (list, tuple, np.ndarray)):
        return [int(x) for x in tag_value]
    if isinstance(tag_value, arr_module.array):
        return [int(x) for x in tag_value]
    if isinstance(tag_value, (bytes, bytearray)):
        n = len(tag_value)
        if n == 0:
            return []
        if n % 4 == 0:
            vals = list(struct.unpack(f'<{n // 4}I', tag_value))
            if all(0 <= v < 100000 for v in vals):
                return vals
        if n % 2 == 0:
            vals = list(struct.unpack(f'<{n // 2}H', tag_value))
            if all(0 <= v < 100000 for v in vals):
                return vals
        if n % 4 == 0:
            return list(struct.unpack(f'<{n // 4}i', tag_value))
        raise ValueError(f'Cannot parse bytes tag of length {n}')
    if isinstance(tag_value, str):
        return [int(x) for x in tag_value.split(',') if x.strip()]
    return [int(tag_value)]


def parse_variant_tag(pv_value):
    if pv_value == "WT":
        return ["WT"]
    try:
        variants = json.loads(pv_value)
        return variants if isinstance(variants, list) else [str(variants)]
    except (json.JSONDecodeError, TypeError):
        return [str(pv_value)]


def parse_variant_id_fields(vid):
    """Parse a variant ID string into position, ref, alt, change_type.

    Handles formats like:
      '6120:t>C'     -> (6120, 't', 'C', 'snv')
      '6120:+T'      -> (6120, '', 'T', 'ins')
      '6120:1c'      -> (6120, 'c', '', 'del')
    Returns (position, ref, alt, change_type) or (None, None, None, None).
    """
    match = re.match(r'(\d+):(.*)', vid)
    if not match:
        return None, None, None, None
    pos = int(match.group(1))
    change = match.group(2)

    # SNV: ref>alt
    snv = re.match(r'([ACGTacgt])>([ACGTacgt])', change)
    if snv:
        return pos, snv.group(1), snv.group(2), 'snv'

    # Insertion: +base(s)
    ins = re.match(r'\+([ACGTacgt]+)', change)
    if ins:
        return pos, '', ins.group(1), 'ins'

    # Deletion: Nbase (e.g., 1c = 1bp deletion of c)
    dele = re.match(r'(\d+)([ACGTacgt])', change)
    if dele:
        return pos, dele.group(2), '', 'del'

    return pos, None, None, 'unknown'


# ============================================================================
# FiberHMM MA-tag parsing
# ============================================================================

def parse_ma_tag(ma_str):
    """Parse the FiberHMM MA:Z structured footprint string.

    Format: "<qlen>;<label>:<s-l>,<s-l>,...;<label>:...;..." where each
    label is a track name optionally followed by a quality-encoding
    suffix (e.g. 'nuc+Q', 'msp+', 'tf+QQQ') and each segment is
    'start-length' with a 1-based query start.

    Returns (qlen, tracks) where tracks maps track name ->
    (starts_0based int32, lengths int32). Starts are converted to
    0-based to match the legacy `ns` convention.
    Robust to empty/missing/malformed sections.
    """
    tracks = {}
    if not ma_str:
        return None, tracks
    parts = ma_str.split(';')
    try:
        qlen = int(parts[0])
    except (ValueError, IndexError):
        qlen = None
    for seg in parts[1:]:
        label, sep, body = seg.partition(':')
        if not sep:
            continue
        track = label.split('+', 1)[0].strip()
        if not track:
            continue
        starts = []
        lengths = []
        if body:
            for pair in body.split(','):
                d = pair.find('-')
                if d <= 0:
                    continue
                try:
                    s = int(pair[:d])
                    ln = int(pair[d + 1:])
                except ValueError:
                    continue
                starts.append(s - 1)  # 1-based -> 0-based query coord
                lengths.append(ln)
        tracks[track] = (np.asarray(starts, dtype=np.int32),
                         np.asarray(lengths, dtype=np.int32))
    return qlen, tracks


def get_read_tracks(read):
    """Per-read footprint tracks as {track: (starts_0based, lengths)}.

    Prefers the FiberHMM MA:Z tag (authoritative; only source of the
    `tf` track). Falls back to the legacy fibertools binary tags when
    MA is absent: ns/nl -> 'nuc', as/al -> 'msp' (no legacy 'tf').
    Always returns int32 numpy arrays; missing tracks are absent.
    """
    try:
        ma = read.get_tag(TAG_MA)
    except KeyError:
        ma = None
    if ma:
        _, tracks = parse_ma_tag(ma)
        return tracks
    # Legacy fallback
    tracks = {}
    for track, (st_tag, ln_tag) in LEGACY_TRACK_TAGS.items():
        try:
            st = read.get_tag(st_tag)
            ln = read.get_tag(ln_tag)
        except KeyError:
            continue
        tracks[track] = (np.asarray(parse_tag_array(st), dtype=np.int32),
                         np.asarray(parse_tag_array(ln), dtype=np.int32))
    return tracks


# ============================================================================
# ReadData
# ============================================================================

class ReadData:
    """Stores per-read coverage matrices and metadata for analysis.

    Supports two construction modes:
    - Pre-allocated (new): pass analysis_length and estimated_reads for
      in-place writes via add_read_coverage(). ~15x less memory when
      analysis_length << ref_length.
    - Legacy (append): pass ref_length only, use add_read() which appends
      rows to lists and vstacks on freeze().
    """

    def __init__(self, ref_length, bin_labels, analysis_length=None,
                 analysis_start=0, estimated_reads=100_000):
        self.ref_length = ref_length
        self.bin_labels = list(bin_labels)
        self.n_bins = len(bin_labels)
        self.bin_to_idx = {label: i for i, label in enumerate(bin_labels)}
        self.n_reads = 0
        self._frozen = False

        # Analysis region metadata
        self.analysis_start = analysis_start
        self.analysis_length = analysis_length if analysis_length else ref_length

        # Pre-allocated mode: matrices sized to analysis_length
        if analysis_length is not None:
            self._preallocated = True
            self._capacity = estimated_reads
            self._nuc_counts = np.empty(self._capacity, dtype=np.int32)
            self._variant_ids = np.empty(self._capacity, dtype=object)
            self._variant_counts = np.empty(self._capacity, dtype=np.int32)
            self._raw_calls = np.empty(self._capacity, dtype=object)
            self._coverage_matrices = {
                label: np.zeros((self._capacity, self.analysis_length),
                                dtype=np.uint8)
                for label in bin_labels
            }
        else:
            # Legacy append mode
            self._preallocated = False
            self._nuc_counts = []
            self._variant_ids = []
            self._variant_counts = []
            self._raw_calls = []
            self._coverage_rows = {label: [] for label in bin_labels}

        # Populated on freeze()
        self.nuc_counts = None
        self.variant_ids = None
        self.variant_counts = None
        self.raw_calls = None   # per-read tuple of raw variant IDs (PR)
        self.coverage_matrices = None
        self.wt_indices = None
        self.variant_indices = None

    def _grow(self):
        """Double capacity of pre-allocated arrays."""
        new_cap = self._capacity * 2
        logging.debug(f"ReadData: growing capacity {self._capacity:,} -> "
                      f"{new_cap:,}")
        new_nc = np.empty(new_cap, dtype=np.int32)
        new_nc[:self._capacity] = self._nuc_counts
        self._nuc_counts = new_nc

        new_vi = np.empty(new_cap, dtype=object)
        new_vi[:self._capacity] = self._variant_ids
        self._variant_ids = new_vi

        new_vc = np.empty(new_cap, dtype=np.int32)
        new_vc[:self._capacity] = self._variant_counts
        self._variant_counts = new_vc

        new_rc = np.empty(new_cap, dtype=object)
        new_rc[:self._capacity] = self._raw_calls
        self._raw_calls = new_rc

        for label in self.bin_labels:
            new_mat = np.zeros((new_cap, self.analysis_length), dtype=np.uint8)
            new_mat[:self._capacity] = self._coverage_matrices[label]
            self._coverage_matrices[label] = new_mat
        self._capacity = new_cap

    def add_read_coverage(self, nuc_count, variant_id, variant_count,
                          bin_coverage, raw_calls=()):
        """Add a read using pre-computed per-bin coverage arrays.

        Parameters
        ----------
        nuc_count : int or None
        variant_id : str
        variant_count : int
        bin_coverage : dict of {bin_label: np.ndarray[uint8]}
            Coverage arrays of length analysis_length, as returned by
            convert_footprints_to_coverage(). Missing bins are all-zeros.
        """
        assert self._preallocated and not self._frozen
        if self.n_reads >= self._capacity:
            self._grow()
        i = self.n_reads
        self._nuc_counts[i] = nuc_count if nuc_count is not None else -1
        self._variant_ids[i] = variant_id
        self._variant_counts[i] = variant_count
        self._raw_calls[i] = raw_calls
        for label, cov in bin_coverage.items():
            if label in self.bin_to_idx:
                self._coverage_matrices[label][i] = cov
        # Bins not in bin_coverage stay as zeros (initialized in __init__)
        self.n_reads += 1

    def add_read(self, nuc_count, variant_id, variant_count, footprints,
                 ref_length, raw_calls=()):
        """Legacy add_read: append-based, used when analysis_length not set."""
        assert not self._preallocated and not self._frozen
        self._nuc_counts.append(nuc_count if nuc_count is not None else -1)
        self._variant_ids.append(variant_id)
        self._variant_counts.append(variant_count)
        self._raw_calls.append(raw_calls)
        coverage = {label: np.zeros(ref_length, dtype=np.uint8)
                    for label in self.bin_labels}
        for bin_label, ref_positions, quality, size in footprints:
            if bin_label not in self.bin_to_idx:
                continue
            pos_arr = np.array(list(ref_positions), dtype=np.int32)
            if len(pos_arr) > 0:
                coverage[bin_label][pos_arr] = 1
        for label in self.bin_labels:
            self._coverage_rows[label].append(coverage[label])
        self.n_reads += 1

    def freeze(self):
        if self._frozen:
            return
        logging.info(f"Freezing ReadData: {self.n_reads:,} reads x "
                     f"{self.analysis_length} bp x {self.n_bins} bins")
        t0 = time.time()

        if self._preallocated:
            n = self.n_reads
            # Trim to actual size with .copy() to release excess memory
            self.nuc_counts = self._nuc_counts[:n].copy()
            self.variant_ids = self._variant_ids[:n].copy()
            self.variant_counts = self._variant_counts[:n].copy()
            self.raw_calls = self._raw_calls[:n].copy()
            self.coverage_matrices = {}
            for label in self.bin_labels:
                self.coverage_matrices[label] = (
                    self._coverage_matrices[label][:n].copy())
                self._coverage_matrices[label] = None
            self._coverage_matrices = None
        else:
            self.nuc_counts = np.array(self._nuc_counts, dtype=np.int32)
            self.variant_ids = np.array(self._variant_ids, dtype=object)
            self.variant_counts = np.array(self._variant_counts, dtype=np.int32)
            self.raw_calls = np.empty(len(self._raw_calls), dtype=object)
            self.raw_calls[:] = self._raw_calls
            self.coverage_matrices = {}
            for label in self.bin_labels:
                if self._coverage_rows[label]:
                    self.coverage_matrices[label] = np.vstack(
                        self._coverage_rows[label])
                else:
                    self.coverage_matrices[label] = np.zeros(
                        (0, self.ref_length), dtype=np.uint8)
                self._coverage_rows[label] = None

        self.wt_indices = np.where(self.variant_ids == 'WT')[0]
        unique_variants = set(self.variant_ids) - {'WT'}
        self.variant_indices = {
            vid: np.where(self.variant_ids == vid)[0]
            for vid in unique_variants
        }
        total_bytes = sum(m.nbytes for m in self.coverage_matrices.values())
        logging.info(f"  Coverage matrices: {total_bytes / 1e6:.1f} MB")
        logging.info(f"  WT reads: {len(self.wt_indices):,}")
        logging.info(f"  Unique variant IDs: {len(self.variant_indices):,}")
        logging.info(f"  Freeze time: {time.time() - t0:.2f}s")
        self._frozen = True
        self._nuc_counts = None
        self._variant_ids = None
        self._variant_counts = None
        self._raw_calls = None

    def get_indices_filtered(self, indices, min_nuc=None, max_nuc=None):
        if min_nuc is None and max_nuc is None:
            return indices
        mask = np.ones(len(indices), dtype=bool)
        nc = self.nuc_counts[indices]
        if min_nuc is not None:
            mask &= (nc >= 0) & (nc >= min_nuc)
        if max_nuc is not None:
            mask &= (nc >= 0) & (nc <= max_nuc)
        return indices[mask]

    def occupancy(self, indices, bin_label):
        if len(indices) == 0:
            return np.zeros(self.analysis_length, dtype=np.float64)
        return self.coverage_matrices[bin_label][indices].mean(axis=0)


# ============================================================================
# BAM Parsing
# ============================================================================

def get_bam_ref_info(bam_path):
    """Read reference sequence info from BAM header.

    Returns dict of {ref_name: ref_length}.
    """
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        return {sq['SN']: sq['LN'] for sq in bam.header['SQ']}


def resolve_target_chrom(ref_info, target_chrom):
    if target_chrom is not None:
        if target_chrom not in ref_info:
            raise ValueError(f"'{target_chrom}' not in BAM. "
                             f"Available: {', '.join(ref_info)}")
        return target_chrom, ref_info[target_chrom]
    if len(ref_info) == 1:
        name = list(ref_info.keys())[0]
        return name, ref_info[name]
    raise ValueError(f"Multiple refs: {', '.join(ref_info)}. "
                     f"Use --target-region.")


def build_query_to_ref_map(aligned_pairs):
    return {qpos: rpos for qpos, rpos in aligned_pairs
            if qpos is not None and rpos is not None}


def convert_footprints(fp_starts, fp_lengths, fp_quals, q2r_map,
                       ref_length, size_bins):
    footprints = []
    for i, (fp_start, fp_len) in enumerate(zip(fp_starts, fp_lengths)):
        bin_label = None
        for min_s, max_s, label in size_bins:
            if min_s <= fp_len <= max_s:
                bin_label = label
                break
        if bin_label is None:
            continue
        qual = (fp_quals[i] if fp_quals is not None and i < len(fp_quals)
                else None)
        ref_positions = set()
        for qpos in range(fp_start, fp_start + fp_len):
            rpos = q2r_map.get(qpos)
            if rpos is not None and 0 <= rpos < ref_length:
                ref_positions.add(rpos)
        if not ref_positions:
            continue
        footprints.append((bin_label, ref_positions, qual, fp_len))
    return footprints


# --- Optimized replacements for build_query_to_ref_map / convert_footprints --

# CIGAR operation codes (from BAM spec)
_CIGAR_M = 0   # alignment match (should not appear in pbmm2 output)
_CIGAR_I = 1   # insertion to reference
_CIGAR_D = 2   # deletion from reference
_CIGAR_N = 3   # skipped region
_CIGAR_S = 4   # soft clip
_CIGAR_H = 5   # hard clip
_CIGAR_P = 6   # padding
_CIGAR_EQ = 7  # sequence match (=)
_CIGAR_X = 8   # sequence mismatch (X)

# Sets for fast lookup
_CONSUMES_QUERY = {_CIGAR_M, _CIGAR_I, _CIGAR_S, _CIGAR_EQ, _CIGAR_X}
_CONSUMES_REF = {_CIGAR_M, _CIGAR_D, _CIGAR_N, _CIGAR_EQ, _CIGAR_X}


def build_q2r_array(read):
    """Build a query-to-reference coordinate mapping array via CIGAR walk.

    Returns a numpy int32 array of length query_length where
    q2r[qpos] = rpos for aligned positions, and -1 for unaligned
    (insertions, soft clips).

    ~7x faster than get_aligned_pairs() + dict comprehension.
    """
    query_length = read.query_length
    q2r = np.full(query_length, -1, dtype=np.int32)
    qpos = 0
    rpos = read.reference_start
    for op, length in read.cigartuples:
        if op in _CONSUMES_QUERY and op in _CONSUMES_REF:
            # M, =, X: both query and ref advance together
            end_q = qpos + length
            q2r[qpos:end_q] = np.arange(rpos, rpos + length, dtype=np.int32)
            qpos = end_q
            rpos += length
        elif op in _CONSUMES_QUERY:
            # I, S: query advances, ref does not
            qpos += length
        elif op in _CONSUMES_REF:
            # D, N: ref advances, query does not
            rpos += length
        # H, P: neither advances
    return q2r


def convert_footprints_to_coverage(fp_starts, fp_lengths, q2r_arr,
                                   analysis_start, analysis_length,
                                   bin_boundaries, bin_labels):
    """Convert footprint tags directly into per-bin coverage arrays.

    Uses numpy fancy indexing on the q2r array instead of per-position
    Python dict lookups. Returns {bin_label: np.uint8 array of length
    analysis_length} with 1 at covered positions, 0 elsewhere.

    Parameters
    ----------
    fp_starts : array-like of int
        Footprint start positions in query coordinates (ns tag).
    fp_lengths : array-like of int
        Footprint lengths (nl tag).
    q2r_arr : np.ndarray[int32]
        Query-to-reference map from build_q2r_array().
    analysis_start : int
        Start of the analysis region in reference coordinates.
    analysis_length : int
        Length of the analysis region.
    bin_boundaries : list of (min_size, max_size)
        Size boundaries for each bin, parallel to bin_labels.
    bin_labels : list of str
        Bin label names, parallel to bin_boundaries.

    Returns
    -------
    dict of {str: np.ndarray[uint8]}
        Coverage arrays per bin label, only for bins that have coverage.
        Missing bins should be treated as all-zeros.
    """
    coverage = {}
    analysis_end = analysis_start + analysis_length
    q2r_len = len(q2r_arr)

    for fp_start, fp_len in zip(fp_starts, fp_lengths):
        # Classify into size bin
        bin_label = None
        for (min_s, max_s), label in zip(bin_boundaries, bin_labels):
            if min_s <= fp_len <= max_s:
                bin_label = label
                break
        if bin_label is None:
            continue

        # Vectorized q2r lookup for footprint span
        fp_end = fp_start + fp_len
        # Clip to query bounds (shouldn't be needed, but safety)
        fp_start_c = max(0, fp_start)
        fp_end_c = min(q2r_len, fp_end)
        if fp_start_c >= fp_end_c:
            continue

        rpositions = q2r_arr[fp_start_c:fp_end_c]
        # Filter: aligned (not -1) and within analysis region
        valid = rpositions[rpositions >= 0]
        valid = valid[(valid >= analysis_start) & (valid < analysis_end)]
        if len(valid) == 0:
            continue

        # Write into coverage array (lazily allocated)
        if bin_label not in coverage:
            coverage[bin_label] = np.zeros(analysis_length, dtype=np.uint8)
        coverage[bin_label][valid - analysis_start] = 1

    return coverage


def convert_tracks_to_coverage(tracks, q2r_arr, analysis_start,
                               analysis_length, analysis_bins):
    """Track-aware coverage builder (FiberHMM).

    Parameters
    ----------
    tracks : dict {track: (starts_0based int32, lengths int32)}
        From get_read_tracks(read).
    q2r_arr : np.ndarray[int32]  query->ref map (build_q2r_array).
    analysis_start, analysis_length : int  reference analysis window.
    analysis_bins : list of (label, track, min_len, max_len)
        max_len inclusive (may be np.inf).

    Returns {label: uint8 array of length analysis_length}. Bins whose
    source track is absent or contributes no covered positions are
    omitted (treated as all-zeros downstream).
    """
    coverage = {}
    analysis_end = analysis_start + analysis_length
    q2r_len = len(q2r_arr)

    for label, track, min_len, max_len in analysis_bins:
        seg = tracks.get(track)
        if seg is None:
            continue
        starts, lengths = seg
        for fp_start, fp_len in zip(starts, lengths):
            if not (min_len <= fp_len <= max_len):
                continue
            fp_start_c = max(0, int(fp_start))
            fp_end_c = min(q2r_len, int(fp_start) + int(fp_len))
            if fp_start_c >= fp_end_c:
                continue
            rpositions = q2r_arr[fp_start_c:fp_end_c]
            valid = rpositions[rpositions >= 0]
            valid = valid[(valid >= analysis_start) & (valid < analysis_end)]
            if len(valid) == 0:
                continue
            if label not in coverage:
                coverage[label] = np.zeros(analysis_length, dtype=np.uint8)
            coverage[label][valid - analysis_start] = 1

    return coverage


def parse_bam(bam_path, analysis_bins, target_chrom=None,
              analysis_region=None, estimated_reads=None):
    """Parse a tagged BAM into a ReadData object (track-aware).

    Parameters
    ----------
    bam_path : str
        Path to BAM file (must have PV/VC tags from stage 2).
    analysis_bins : list of (label, track, min_len, max_len)
        Track-aware analysis bins (see build_analysis_bins). Footprint
        segments are read per-track via get_read_tracks() (FiberHMM
        MA:Z, or legacy ns/nl/as fallback).
    target_chrom : str or None
        Target chromosome name (auto-detected if single-ref BAM).
    analysis_region : tuple of (start, end) or None
        Reference window for coverage matrices. None = full reference.
    estimated_reads : int or None
        Hint for pre-allocation size. If None, tries to read from index.

    Returns
    -------
    rd : ReadData
    ref_length : int
    ref_name : str
    stats : dict
    """
    stats = {'total': 0, 'unmapped': 0, 'missing_pv_vc': 0,
             'missing_fp_tags': 0, 'no_footprints': 0, 'off_target': 0,
             'parsed': 0}
    bin_labels = [b[0] for b in analysis_bins]

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        ref_info = {sq['SN']: sq['LN'] for sq in bam.header['SQ']}
        logging.info(f"References in BAM: {ref_info}")
        ref_name, ref_length = resolve_target_chrom(ref_info, target_chrom)
        logging.info(f"Using reference: {ref_name} ({ref_length:,} bp)")

        if estimated_reads is None:
            try:
                estimated_reads = sum(s.mapped for s in
                                      bam.get_index_statistics())
            except (ValueError, AttributeError):
                estimated_reads = 100_000

        total_reads = None
        try:
            total_reads = sum(s.mapped + s.unmapped
                              for s in bam.get_index_statistics())
        except (ValueError, AttributeError):
            pass

        # Single coverage path. Window = analysis_region or full ref.
        if analysis_region is not None:
            a_start, a_end = analysis_region
        else:
            a_start, a_end = 0, ref_length
        a_length = a_end - a_start
        logging.info(f"Analysis region: {a_start}-{a_end} "
                     f"({a_length} bp of {ref_length} bp); "
                     f"bins={bin_labels}")
        rd = ReadData(ref_length, bin_labels,
                      analysis_length=a_length,
                      analysis_start=a_start,
                      estimated_reads=estimated_reads)

        for read in tqdm(bam, total=total_reads, desc="Parsing BAM",
                         unit=" reads", disable=not HAS_TQDM):
            stats['total'] += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                stats['unmapped'] += 1
                continue
            if read.reference_name != ref_name:
                stats['off_target'] += 1
                continue
            try:
                pv_value = read.get_tag(TAG_PROMOTER_VARIANT)
                vc_value = read.get_tag(TAG_VARIANT_COUNT)
            except KeyError:
                stats['missing_pv_vc'] += 1
                continue
            try:
                nuc_count = int(read.get_tag(TAG_NUC_COUNT))
            except KeyError:
                nuc_count = None
            # Per-read RAW (pre-cluster-consensus) calls for the
            # co-occupancy secondary-mutation control (§2.5). Absent on
            # older BAMs -> empty (feature simply has no data then).
            try:
                pr_value = read.get_tag(TAG_RAW_VARIANT)
            except KeyError:
                pr_value = None
            if pr_value is None or pr_value == "WT":
                raw_calls = ()
            else:
                raw_calls = tuple(
                    v for v in parse_variant_tag(pr_value) if v != "WT")

            tracks = get_read_tracks(read)
            if not tracks:
                stats['missing_fp_tags'] += 1
                continue
            if sum(len(s[0]) for s in tracks.values()) == 0:
                stats['no_footprints'] += 1

            if pv_value == "WT" or vc_value == 0:
                variant_id = "WT"
            else:
                variants = parse_variant_tag(pv_value)
                variant_id = (variants[0] if len(variants) == 1
                              else json.dumps(variants))

            q2r_arr = build_q2r_array(read)
            bin_cov = convert_tracks_to_coverage(
                tracks, q2r_arr, a_start, a_length, analysis_bins)
            rd.add_read_coverage(nuc_count, variant_id, int(vc_value),
                                 bin_cov, raw_calls=raw_calls)
            stats['parsed'] += 1

    logging.info(f"Parsed {stats['parsed']:,} / {stats['total']:,}")
    rd.freeze()
    return rd, ref_length, ref_name, stats


# ============================================================================
# Variant Grouping
# ============================================================================

def group_variants(rd, include_multi=False, min_reads=50,
                   min_nuc=None, max_nuc=None):
    """Identify testable variants and their read indices."""
    logging.info("Grouping variants...")

    single_variants = {}
    multi_variants = {}

    for raw_vid, raw_indices in rd.variant_indices.items():
        parsed = parse_variant_tag(raw_vid)
        if len(parsed) == 1 and parsed[0] != 'WT':
            canonical = parsed[0]
            if canonical not in single_variants:
                single_variants[canonical] = []
            single_variants[canonical].append(raw_indices)
        elif len(parsed) > 1:
            multi_variants[raw_vid] = (parsed, raw_indices)

    variant_indices = {}
    for vid, idx_list in single_variants.items():
        variant_indices[vid] = np.concatenate(idx_list)

    multi_assigned = 0
    if include_multi:
        for raw_vid, (parsed, raw_indices) in multi_variants.items():
            primary = parsed[0]
            if primary not in variant_indices:
                variant_indices[primary] = raw_indices
            else:
                variant_indices[primary] = np.concatenate(
                    [variant_indices[primary], raw_indices])
            multi_assigned += len(raw_indices)

    filtered = {}
    skipped_low = 0
    for vid, indices in variant_indices.items():
        f_idx = rd.get_indices_filtered(indices, min_nuc=min_nuc,
                                         max_nuc=max_nuc)
        if len(f_idx) >= min_reads:
            filtered[vid] = f_idx
        else:
            skipped_low += 1

    logging.info(f"  Single-variant IDs: {len(single_variants)}")
    logging.info(f"  Multi-variant raw IDs: {len(multi_variants)}")
    if include_multi:
        logging.info(f"  Multi-variant reads assigned: {multi_assigned:,}")
    logging.info(f"  Variants passing min_reads={min_reads}: {len(filtered)}")
    logging.info(f"  Variants skipped (low reads): {skipped_low}")

    if filtered:
        counts = sorted([len(v) for v in filtered.values()])
        logging.info(f"  Coverage distribution: "
                     f"min={counts[0]}, "
                     f"p25={counts[len(counts)//4]}, "
                     f"median={counts[len(counts)//2]}, "
                     f"p75={counts[3*len(counts)//4]}, "
                     f"max={counts[-1]}")

    return filtered


# ============================================================================
# Ground Truth NC & NC-Matched Subsampling
# ============================================================================

def compute_ground_truth_nc(rd, wt_idx, min_nuc=None, max_nuc=None):
    wt_f = rd.get_indices_filtered(wt_idx, min_nuc=min_nuc, max_nuc=max_nuc)
    nc = rd.nuc_counts[wt_f]
    nc_valid = nc[nc >= 0]
    if len(nc_valid) == 0:
        return None
    unique_nc, counts = np.unique(nc_valid, return_counts=True)
    fracs = counts.astype(np.float64) / counts.sum()
    logging.info(f"  Ground truth NC: {len(wt_f):,} reads, "
                 f"range=[{unique_nc[0]}, {unique_nc[-1]}], "
                 f"mean={np.mean(nc_valid):.1f}")
    return {'nc_vals': unique_nc, 'nc_fracs': fracs, 'n_reads': len(wt_f)}


def nc_matched_subsample(indices, nuc_counts, target_nc_fracs,
                         target_nc_vals, rng, replace=False):
    nc = nuc_counts[indices]
    available = {nc_val: indices[nc == nc_val] for nc_val in target_nc_vals}
    if not replace:
        max_n = len(indices)
        for nc_val, frac in zip(target_nc_vals, target_nc_fracs):
            pool_size = len(available.get(nc_val, []))
            if frac > 1e-10 and pool_size > 0:
                max_n = min(max_n, int(pool_size / frac))
            elif frac > 1e-10:
                max_n = 0
        if max_n == 0:
            return indices, len(indices)
        draws = []
        for nc_val, frac in zip(target_nc_vals, target_nc_fracs):
            pool = available.get(nc_val, np.array([], dtype=int))
            n_take = min(max(1, int(round(frac * max_n))), len(pool))
            if n_take > 0:
                draws.append(rng.choice(pool, size=n_take, replace=False))
        if not draws:
            return indices, len(indices)
        return np.concatenate(draws), sum(len(d) for d in draws)
    else:
        size = len(indices)
        draws = []
        for nc_val, frac in zip(target_nc_vals, target_nc_fracs):
            pool = available.get(nc_val, np.array([], dtype=int))
            if len(pool) == 0:
                continue
            draws.append(rng.choice(pool, size=max(1, int(round(frac * size))),
                                    replace=True))
        if not draws:
            return rng.choice(indices, size=size, replace=True), size
        combined = np.concatenate(draws)
        if len(combined) >= size:
            return rng.permutation(combined)[:size], size
        pad = rng.choice(combined, size=size - len(combined), replace=True)
        return np.concatenate([combined, pad]), size


def _nc_matched_batch_draw(nc_draw_counts, wt_nc_pools, n_samples,
                           subsample_size, rng):
    all_draws = []
    for nc_val, n_draw in nc_draw_counts:
        pool = wt_nc_pools.get(nc_val)
        if pool is None or len(pool) == 0:
            continue
        all_draws.append(rng.choice(pool, size=(n_samples, n_draw),
                                    replace=True))
    if not all_draws:
        fallback = max(wt_nc_pools.values(), key=len)
        return rng.choice(fallback, size=(n_samples, subsample_size),
                          replace=True)
    combined = np.concatenate(all_draws, axis=1)
    if combined.shape[1] >= subsample_size:
        return combined[:, :subsample_size]
    pad_size = subsample_size - combined.shape[1]
    pad = rng.choice(combined.ravel(), size=(n_samples, pad_size),
                     replace=True)
    return np.concatenate([combined, pad], axis=1)


# ============================================================================
# Shared Memory Helpers
# ============================================================================

def create_shared_matrices(wt_mats, bin_labels):
    shm_info = {}
    shm_objects = []
    for label in bin_labels:
        mat = wt_mats[label]
        shm = shared_memory.SharedMemory(create=True, size=mat.nbytes)
        shared_arr = np.ndarray(mat.shape, dtype=mat.dtype, buffer=shm.buf)
        shared_arr[:] = mat[:]
        shm_info[label] = {'name': shm.name, 'shape': mat.shape,
                           'dtype': str(mat.dtype)}
        shm_objects.append(shm)
    return shm_info, shm_objects


def load_shared_matrix(shm_meta):
    shm = shared_memory.SharedMemory(name=shm_meta['name'], create=False)
    arr = np.ndarray(shm_meta['shape'], dtype=np.dtype(shm_meta['dtype']),
                     buffer=shm.buf)
    return arr, shm


def cleanup_shared_memory(shm_objects):
    for shm in shm_objects:
        try:
            shm.close()
            shm.unlink()
        except Exception:
            pass


# ============================================================================
# Cluster Detection
# ============================================================================

def _detect_clusters_core(sig_mask, abs_delta, signed_delta=None,
                          gap_tolerance=2, merge_distance=5, min_width=3):
    """Unified cluster detector. Behavior-preserving merge of the former
    _detect_clusters_adaptive (signed_delta=None) and
    _detect_clusters_from_mask (signed_delta given -> adds
    mean_signed_delta + direction). Identical numerics to both.
    """
    n = len(sig_mask)
    if not np.any(sig_mask):
        return []
    if gap_tolerance > 0:
        # Vectorized gap-bridge (was a per-position Python np.any loop —
        # the profiled hot path: 40000 calls/variant, ~10M np.any).
        # A position is bridged iff some True lies within gap_tolerance
        # to its left AND within gap_tolerance to its right. Window
        # sums via prefix-sum; bitwise-identical to the loop. (sig
        # positions stay True under the OR, exactly as the copy did.)
        g = gap_tolerance
        c = np.empty(n + 1, dtype=np.int64)
        c[0] = 0
        np.cumsum(sig_mask.astype(np.int64), out=c[1:])
        idx = np.arange(n)
        left_any = (c[idx] - c[np.maximum(0, idx - g)]) > 0
        right_any = (c[np.minimum(n, idx + g + 1)]
                     - c[np.minimum(n, idx + 1)]) > 0
        sig_mask = sig_mask | (left_any & right_any)
    labeled, n_clusters = ndimage.label(sig_mask)
    raw = []
    for c in range(1, n_clusters + 1):
        pos = np.where(labeled == c)[0]
        raw.append({'start': int(pos[0]), 'end': int(pos[-1])})
    if len(raw) > 1:
        merged = [raw[0].copy()]
        for c in raw[1:]:
            if c['start'] - merged[-1]['end'] <= merge_distance:
                merged[-1]['end'] = c['end']
            else:
                merged.append(c.copy())
        raw = merged
    clusters = []
    for c in raw:
        s, e = c['start'], c['end']
        w = e - s + 1
        if w < min_width:
            continue
        ra = abs_delta[s:e + 1]
        cl = {
            'start': s, 'end': e, 'width': w,
            'sum_abs_delta': float(np.sum(ra)),
            'max_abs_delta': float(np.max(ra)),
            'mean_abs_delta': float(np.mean(ra)),
            'peak_position': s + int(np.argmax(ra)),
        }
        if signed_delta is not None:
            rs = signed_delta[s:e + 1]
            cl['mean_signed_delta'] = float(np.mean(rs))
            cl['direction'] = 'loss' if np.mean(rs) < 0 else 'gain'
        clusters.append(cl)
    return clusters


def _detect_clusters_adaptive(abs_delta, threshold, gap_tolerance=2,
                              merge_distance=5, min_width=3):
    return _detect_clusters_core(abs_delta > threshold, abs_delta, None,
                                 gap_tolerance, merge_distance, min_width)


def _detect_clusters_from_mask(sig_mask, abs_delta, signed_delta,
                               gap_tolerance=2, merge_distance=5,
                               min_width=3):
    return _detect_clusters_core(sig_mask, abs_delta, signed_delta,
                                 gap_tolerance, merge_distance, min_width)


# ============================================================================
# Null Calibration — Worker
# ============================================================================

def _pool_warmup(_):
    """Trivial task to force a ProcessPoolExecutor worker to spawn. The
    short sleep keeps a fast worker from grabbing several tasks before
    its peers start, so mapping n_workers of these spins up all of
    them. Module-level so it pickles."""
    import time as _t
    _t.sleep(0.05)
    return None


def _null_worker_shared(args):
    (it_seeds, shm_info, bin_labels, nc_draw_counts,
     wt_nc_pools_ser, subsample_size, wt_occ_list,
     analysis_length, cluster_threshold_quantile, gap_tolerance,
     merge_distance, absolute_delta_threshold) = args

    n_iters = len(it_seeds)
    n_bins = len(bin_labels)

    shared_mats = {}
    shm_handles = []
    for label in bin_labels:
        arr, shm = load_shared_matrix(shm_info[label])
        shared_mats[label] = arr
        shm_handles.append(shm)

    wt_nc_pools = {v: p for v, p in wt_nc_pools_ser}
    wt_occ = {l: np.array(v) for l, v in zip(bin_labels, wt_occ_list)}

    delta_out = np.zeros((n_iters, n_bins, analysis_length), dtype=np.float32)

    # Workers fill null_delta ONLY. Per-iteration cluster detection used
    # to be computed here too, but the caller discarded it (`dc, _ =`)
    # and recomputed clusters post-hoc with the unified rule. Computing
    # + pickling it was pure waste (profiled hot path). Numerically a
    # no-op: the returned delta_out is unchanged.
    for li, seed in enumerate(it_seeds):
        rng = np.random.default_rng(seed)
        pv_idx = _nc_matched_batch_draw(
            nc_draw_counts, wt_nc_pools, 1, subsample_size, rng)[0]
        for bi, label in enumerate(bin_labels):
            pv_occ = shared_mats[label][pv_idx].mean(axis=0)
            delta_out[li, bi, :] = (
                pv_occ - wt_occ[label]).astype(np.float32)

    for shm in shm_handles:
        shm.close()
    return delta_out


# ============================================================================
# Null Calibration — Main
# ============================================================================

def build_wt_null_context(rd, wt_idx, min_nuc, max_nuc, analysis_region,
                           make_shared=True):
    """Build the WT-side structures that are INVARIANT across variants
    (the WT read set, its coverage matrices sliced to the analysis
    window, NC pools, and optionally a shared-memory copy of the
    matrices). These depend only on (rd, wt_idx, filters,
    analysis_region) — not on any variant — so for the per-variant null
    they can be built ONCE and reused across all variants instead of
    rebuilt 149x. Extracted verbatim from run_null_calibration; the
    arrays produced are identical to the per-call build.
    """
    bin_labels = rd.bin_labels
    if analysis_region is not None:
        a_start, a_end = analysis_region
    else:
        a_start, a_end = 0, rd.ref_length
    analysis_length = a_end - a_start

    wt_f = rd.get_indices_filtered(wt_idx, min_nuc=min_nuc, max_nuc=max_nuc)
    needs_col_slice = (rd.analysis_length != analysis_length
                       or rd.analysis_start != a_start)
    if needs_col_slice:
        col_indices = np.arange(a_start, a_end)
        wt_mats = {l: rd.coverage_matrices[l][np.ix_(wt_f, col_indices)]
                   for l in bin_labels}
    else:
        wt_mats = {l: rd.coverage_matrices[l][wt_f]
                   for l in bin_labels}
    wt_nc = rd.nuc_counts[wt_f]
    wt_nc_pools = {}
    for nc_val in np.unique(wt_nc):
        pool = np.where(wt_nc == nc_val)[0]
        if len(pool) > 0:
            wt_nc_pools[int(nc_val)] = pool
    wt_nc_pools_ser = list(wt_nc_pools.items())

    shm_info, shm_objects = (create_shared_matrices(wt_mats, bin_labels)
                             if make_shared else (None, []))
    return {
        'wt_f': wt_f, 'wt_mats': wt_mats, 'wt_nc': wt_nc,
        'wt_nc_pools': wt_nc_pools, 'wt_nc_pools_ser': wt_nc_pools_ser,
        'shm_info': shm_info, 'shm_objects': shm_objects,
        'analysis_region': (a_start, a_end), 'bin_labels': list(bin_labels),
    }


def run_null_calibration(rd, wt_idx, ground_truth_nc, coverage_level,
                         min_nuc=None, max_nuc=None,
                         n_iterations=2000, random_seed=42,
                         analysis_region=None,
                         cluster_threshold_quantile=0.95,
                         absolute_delta_threshold=None,
                         gap_tolerance=2, merge_distance=5,
                         n_workers=1, wt_ctx=None, executor=None):
    """Generate empirical null at a given coverage level.

    wt_ctx (optional): a build_wt_null_context() dict. When given, the
    invariant WT slice / NC pools / shared memory are reused instead of
    rebuilt — numerically identical (same arrays), just not recomputed
    per variant. executor (optional): a persistent ProcessPoolExecutor
    to reuse instead of creating/tearing down one per call. Neither
    changes any seed, the chunking, the Option-A reference, or the
    null math — only what is rebuilt vs reused.
    """
    logging.info(f"  Null calibration: N={coverage_level}, "
                 f"iters={n_iterations}, workers={n_workers}")

    bin_labels = rd.bin_labels
    n_bins = len(bin_labels)
    if analysis_region is not None:
        a_start, a_end = analysis_region
    else:
        a_start, a_end = 0, rd.ref_length
    analysis_length = a_end - a_start

    # Invariant WT structures: reuse the caller's prebuilt context when
    # given (built once, identical arrays), else build locally — same
    # arrays either way; only whether they are recomputed per variant.
    if wt_ctx is None:
        wt_ctx = build_wt_null_context(rd, wt_idx, min_nuc, max_nuc,
                                       analysis_region, make_shared=False)
    wt_f = wt_ctx['wt_f']
    wt_mats = wt_ctx['wt_mats']
    wt_nc = wt_ctx['wt_nc']
    wt_nc_pools = wt_ctx['wt_nc_pools']
    wt_nc_pools_ser = wt_ctx['wt_nc_pools_ser']
    nc_vals = ground_truth_nc['nc_vals']
    nc_fracs = ground_truth_nc['nc_fracs']

    # Option A reference: WT mean occupancy reweighted to the passed NC
    # distribution (the variant's, for per-variant nulls). This removes
    # the nucleosome-count confounder — both the observed delta and the
    # null delta are referenced to "what WT does at this NC profile",
    # so the null is centered at ~0 under H0. Renormalize over NC values
    # that actually have WT reads.
    wt_occ = {}
    for l in bin_labels:
        acc = np.zeros(analysis_length, dtype=np.float64)
        tot = 0.0
        for nc_val, frac in zip(nc_vals, nc_fracs):
            pool = wt_nc_pools.get(nc_val)
            if pool is None or len(pool) == 0 or frac <= 0:
                continue
            acc += frac * wt_mats[l][pool].mean(axis=0)
            tot += frac
        if tot > 0:
            wt_occ[l] = acc / tot
        else:
            wt_occ[l] = wt_mats[l].mean(axis=0).astype(np.float64)

    nc_draw_counts = []
    for nc_val, frac in zip(nc_vals, nc_fracs):
        n_draw = max(1, int(round(frac * coverage_level)))
        if nc_val in wt_nc_pools:
            nc_draw_counts.append((nc_val, n_draw))

    master_rng = np.random.default_rng(random_seed)
    all_seeds = master_rng.integers(0, 2**63, size=n_iterations)
    wt_occ_list = [wt_occ[l].tolist() for l in bin_labels]

    null_delta = np.zeros((n_iterations, n_bins, analysis_length),
                          dtype=np.float32)
    # ---- Fill null_delta only (cluster detection done post-hoc with the
    # SAME rule used for the observed data — fixes Issue 2b) ----
    if n_workers > 1 and n_iterations > 1:
        # Reuse caller-provided shared memory / executor when present
        # (built once, reused across all variants); else own them for
        # this call. Chunking and seeds are unchanged either way.
        if wt_ctx.get('shm_info') is not None:
            shm_info, shm_objects, _own_shm = (
                wt_ctx['shm_info'], [], False)
        else:
            shm_info, shm_objects = create_shared_matrices(
                wt_mats, bin_labels)
            _own_shm = True
        chunk_size = max(1, n_iterations // n_workers)
        seed_chunks = [all_seeds[i:i + chunk_size]
                       for i in range(0, n_iterations, chunk_size)]
        worker_args = [
            (chunk, shm_info, bin_labels, nc_draw_counts,
             wt_nc_pools_ser, coverage_level, wt_occ_list,
             analysis_length, cluster_threshold_quantile, gap_tolerance,
             merge_distance, absolute_delta_threshold)
            for chunk in seed_chunks
        ]
        try:
            if executor is not None:
                ex, _own_ex = executor, False
            else:
                ex, _own_ex = ProcessPoolExecutor(max_workers=n_workers), True
            try:
                futures = {ex.submit(_null_worker_shared, wa): ci
                           for ci, wa in enumerate(worker_args)}
                results_by_chunk = {}
                for future in as_completed(futures):
                    results_by_chunk[futures[future]] = future.result()
                offset = 0
                for ci in range(len(seed_chunks)):
                    dc = results_by_chunk[ci]
                    n_c = dc.shape[0]
                    null_delta[offset:offset + n_c] = dc
                    offset += n_c
            finally:
                if _own_ex:
                    ex.shutdown()
        finally:
            if _own_shm:
                cleanup_shared_memory(shm_objects)
    else:
        for it in tqdm(range(n_iterations), desc="    Null iters",
                       disable=not HAS_TQDM):
            rng = np.random.default_rng(all_seeds[it])
            pv_idx = _nc_matched_batch_draw(
                nc_draw_counts, wt_nc_pools, 1, coverage_level, rng)[0]
            for bi, label in enumerate(bin_labels):
                pv_occ = wt_mats[label][pv_idx].mean(axis=0)
                null_delta[it, bi, :] = (
                    pv_occ - wt_occ[label]).astype(np.float32)

    # ---- Unified cluster rule + family-wide max-statistic null ----
    results = {
        'null_delta': null_delta, 'wt_occ': wt_occ,
        'coverage_level': coverage_level, 'n_iterations': n_iterations,
        'analysis_region': (a_start, a_end),
        'analysis_length': analysis_length,
        'summary': {}, 'null_cluster_sums': {},
        'null_max_cluster_sums': {}, 'pos_thresh': {},
    }
    abs_null = {}
    for bi, label in enumerate(bin_labels):
        d = null_delta[:, bi, :]
        ad = np.abs(d)
        abs_null[label] = ad
        results['summary'][label] = {
            'pos_null_mean': d.mean(axis=0),
            'pos_null_std': d.std(axis=0),
        }
        if cluster_threshold_quantile is not None:
            pt = np.percentile(ad, cluster_threshold_quantile * 100, axis=0)
        else:
            pt = np.zeros(analysis_length)
        if absolute_delta_threshold is not None:
            pt = np.maximum(pt, absolute_delta_threshold)
        results['pos_thresh'][label] = pt

    # Per iteration: detect clusters in every bin with the per-position
    # threshold above (identical to the observed-data rule), record
    # per-bin sums and the family-wide (across-bin) max cluster Σ|Δ|.
    null_familywise_max = np.zeros(n_iterations, dtype=np.float32)
    per_bin_sums = {l: [] for l in bin_labels}
    per_bin_iter_max = {l: np.zeros(n_iterations, dtype=np.float32)
                        for l in bin_labels}
    for it in range(n_iterations):
        fam = 0.0
        for bi, label in enumerate(bin_labels):
            ad_it = abs_null[label][it]
            sig = ad_it > results['pos_thresh'][label]
            cl = _detect_clusters_from_mask(
                sig, ad_it, null_delta[it, bi, :],
                gap_tolerance, merge_distance)
            csums = [c['sum_abs_delta'] for c in cl]
            per_bin_sums[label].extend(csums)
            m = max(csums, default=0.0)
            per_bin_iter_max[label][it] = m
            if m > fam:
                fam = m
        null_familywise_max[it] = fam
    for label in bin_labels:
        results['null_cluster_sums'][label] = np.asarray(
            per_bin_sums[label], dtype=np.float32)
        results['null_max_cluster_sums'][label] = per_bin_iter_max[label]
    results['null_familywise_max'] = null_familywise_max
    return results


# ============================================================================
# Auto Coverage Grid
# ============================================================================

def auto_coverage_grid(variant_groups, n_points=8):
    """Compute coverage grid from the actual variant read count distribution."""
    counts = sorted([len(v) for v in variant_groups.values()])
    if len(counts) == 0:
        return [50, 100, 200, 500, 1000, 2000]
    percentiles = [5, 10, 25, 50, 75, 90, 95]
    grid = set()
    for p in percentiles:
        idx = int(len(counts) * p / 100)
        idx = min(idx, len(counts) - 1)
        grid.add(counts[idx])
    grid.add(counts[0])
    grid.add(counts[-1])
    nice_grid = sorted(set(max(20, int(round(g / 10) * 10)) for g in grid))
    logging.info(f"  Auto coverage grid ({len(nice_grid)} depths): {nice_grid}")
    return nice_grid


# ============================================================================
# Per-Variant Testing
# ============================================================================

def _compute_variant_result(vid, var_indices, var_matched, var_n,
                            var_occ_fn, null_res, closest, bin_labels,
                            a_start, analysis_length,
                            prom_rel_s, prom_rel_e, prom_slice,
                            cluster_threshold_quantile,
                            absolute_delta_threshold,
                            gap_tolerance, merge_distance,
                            mde_alpha=0.05):
    """Shared per-variant testing core. Verbatim extraction of the
    per-bin block that was duplicated between test_single_variant and
    _variant_test_worker. The only thing that differed between callers
    was how var_occ is obtained, injected here as var_occ_fn(label).
    Numerics are identical to the former duplicated code.
    """
    wt_occ = null_res['wt_occ']
    n_null = null_res['null_delta'].shape[0]
    # Minimum detectable effect: smallest |Δ| this variant's N could
    # resolve at each position given the WT (Option-A) occupancy p.
    # Lets a true null be distinguished from "underpowered".
    z_mde = float(norm.ppf(1.0 - mde_alpha / 2.0))
    mde_all = []
    # Family-wide max-statistic null (computed in run_null_calibration by
    # the SAME detection rule as below): per null iteration, the max
    # cluster Σ|Δ| across all bins. Calibrating the observed family-wide
    # max against this is the correct multiplicity-aware test (fixes
    # Issue 1: no uncorrected min-over-bins; Issue 2: identical null/
    # observed rule, no selection-biased pooled reference).
    null_fwm = null_res.get('null_familywise_max',
                            np.zeros(n_null, dtype=np.float32))
    pos_thresh_by_label = null_res.get('pos_thresh', {})

    results = {'variant_id': vid, 'n_raw': len(var_indices),
               'n_nc_matched': var_n, 'null_depth_used': closest}
    obs_familywise_max = 0.0

    for bi, label in enumerate(bin_labels):
        var_occ = var_occ_fn(label)
        delta_obs = var_occ - wt_occ[label]

        null_delta = null_res['null_delta'][:, bi, :]
        abs_obs = np.abs(delta_obs)
        abs_null = np.abs(null_delta)

        # Per-position empirical p / BH-q — DIAGNOSTIC ONLY (positions
        # are strongly spatially autocorrelated; not a calibrated count).
        exceed = np.sum(abs_null >= abs_obs[np.newaxis, :], axis=0)
        empirical_p = (exceed + 1) / (n_null + 1)
        prom_p = empirical_p[prom_slice]
        prom_q = benjamini_hochberg(prom_p)
        q_values = np.ones(analysis_length)
        q_values[prom_slice] = prom_q

        null_mean = null_res['summary'][label]['pos_null_mean']
        null_std = np.maximum(null_res['summary'][label]['pos_null_std'], 1e-6)
        z_scores = (delta_obs - null_mean) / null_std

        # Observed clusters via the EXACT threshold the null used.
        pos_thresh = pos_thresh_by_label.get(label)
        if pos_thresh is None:
            pos_thresh = np.percentile(
                abs_null, cluster_threshold_quantile * 100, axis=0)
        sig_mask = abs_obs > pos_thresh
        if absolute_delta_threshold is not None:
            sig_mask &= abs_obs > absolute_delta_threshold

        clusters = _detect_clusters_from_mask(
            sig_mask, abs_obs, delta_obs, gap_tolerance, merge_distance)

        # Bin contribution to the observed family-wide max (mirrors the
        # null computation: max cluster Σ|Δ| over all detected clusters).
        for c in clusters:
            if c['sum_abs_delta'] > obs_familywise_max:
                obs_familywise_max = c['sum_abs_delta']

        null_max_csums = null_res['null_max_cluster_sums'].get(
            label, np.array([]))
        sig_clusters = []
        for c in clusters:
            if c['end'] < prom_rel_s or c['start'] >= prom_rel_e:
                continue
            # Per-bin diagnostic p (per-iteration max within this bin).
            if len(null_max_csums) > 0:
                c['max_sum_p'] = max(
                    float(np.mean(null_max_csums >= c['sum_abs_delta'])),
                    1.0 / (len(null_max_csums) + 1))
            else:
                c['max_sum_p'] = np.nan
            c['abs_start'] = c['start'] + a_start
            c['abs_end'] = c['end'] + a_start
            if c.get('max_sum_p', 1.0) < 0.05:
                sig_clusters.append(c)

        n_sig_pos = int(np.sum(q_values[prom_slice] < 0.10))

        p_wt = np.clip(wt_occ[label], 0.0, 1.0)
        mde = z_mde * np.sqrt(p_wt * (1.0 - p_wt) / max(var_n, 1))
        mde_all.append(mde[prom_slice])

        results[label] = {
            'delta_obs': delta_obs.astype(np.float32),
            'variant_occ': var_occ.astype(np.float32),
            'empirical_p': empirical_p.astype(np.float32),
            'q_values': q_values.astype(np.float32),
            'z_scores': z_scores.astype(np.float32),
            'mde': mde.astype(np.float32),
            'n_sig_positions_fdr10': n_sig_pos,
            'max_abs_delta': float(np.max(abs_obs)),
            'significant_clusters': sig_clusters,
            'all_promoter_clusters': [c for c in clusters
                                      if c['end'] >= prom_rel_s
                                      and c['start'] < prom_rel_e],
        }

    # Calibrated, multiplicity-aware variant p-value: where does the
    # observed family-wide max fall in the null family-wide-max dist?
    results['obs_familywise_max'] = float(obs_familywise_max)
    results['best_cluster_p'] = float(
        (np.sum(null_fwm >= obs_familywise_max) + 1) / (n_null + 1))
    results['mde_median'] = (float(np.median(np.concatenate(mde_all)))
                             if mde_all else float('nan'))
    return results


def test_single_variant(rd, wt_idx, var_indices, var_id,
                         analysis_region, promoter_start, promoter_end,
                         min_nuc=None, max_nuc=None, random_seed=42,
                         n_null_iterations=2000,
                         cluster_threshold_quantile=0.95,
                         absolute_delta_threshold=None,
                         gap_tolerance=2, merge_distance=5):
    """Test one variant against its OWN per-variant null (Issue 3 fix).

    The null is built at the variant's exact read count N, sampling WT
    WITH replacement, NC-matched to the variant's OWN NC distribution,
    and referenced to the Option-A NC-reweighted WT mean (computed
    inside run_null_calibration). The FULL variant read set is used as
    a point estimate (no NC-subsampling of the variant), so there is no
    null/observed sampling-scheme mismatch.
    """
    bin_labels = rd.bin_labels
    a_start, a_end = analysis_region
    analysis_length = a_end - a_start
    prom_rel_s = max(0, promoter_start - a_start)
    prom_rel_e = min(analysis_length, promoter_end - a_start)
    prom_slice = slice(prom_rel_s, prom_rel_e)

    var_f = rd.get_indices_filtered(var_indices, min_nuc=min_nuc,
                                    max_nuc=max_nuc)
    var_n = len(var_f)
    if var_n < 20:
        return None

    # The variant's own NC distribution drives both the null NC-matching
    # and the Option-A reference.
    variant_nc = compute_ground_truth_nc(rd, var_f, min_nuc=min_nuc,
                                         max_nuc=max_nuc)
    if variant_nc is None:
        return None

    null_res = run_null_calibration(
        rd, wt_idx, variant_nc, coverage_level=var_n,
        min_nuc=min_nuc, max_nuc=max_nuc,
        n_iterations=n_null_iterations, random_seed=random_seed,
        analysis_region=analysis_region,
        cluster_threshold_quantile=cluster_threshold_quantile,
        absolute_delta_threshold=absolute_delta_threshold,
        gap_tolerance=gap_tolerance, merge_distance=merge_distance,
        n_workers=1)

    needs_col_slice = (rd.analysis_length != analysis_length
                       or rd.analysis_start != a_start)
    if needs_col_slice:
        col_indices = np.arange(a_start, a_end)

        def var_occ_fn(label):
            return rd.coverage_matrices[label][
                np.ix_(var_f, col_indices)].mean(axis=0)
    else:
        def var_occ_fn(label):
            return rd.coverage_matrices[label][var_f].mean(axis=0)

    return _compute_variant_result(
        var_id, var_indices, var_f, var_n, var_occ_fn, null_res,
        var_n, bin_labels, a_start, analysis_length,
        prom_rel_s, prom_rel_e, prom_slice,
        cluster_threshold_quantile, absolute_delta_threshold,
        gap_tolerance, merge_distance)


# ============================================================================
# Per-Variant Testing (per-variant nulls; Issue 3 fix)
# ============================================================================

def _nc_wasserstein(nc_a, nc_b):
    """Wasserstein-1 between two NC distributions (compute_ground_truth_nc
    dicts: nc_vals as support, nc_fracs as weights)."""
    return float(wasserstein_distance(
        nc_a['nc_vals'], nc_b['nc_vals'],
        nc_a['nc_fracs'], nc_b['nc_fracs']))


def _w1_nonneg_int(a, b):
    """Exact 1D Wasserstein-1 between two empirical samples of
    NON-NEGATIVE INTEGERS. W1 = ∫|F_a − F_b| dx; for integer support
    F is constant on each unit interval, so this equals
    Σ_k |F_a(k) − F_b(k)| over k = 0..max. This is mathematically the
    same value scipy.stats.wasserstein_distance returns for integer
    samples, but O(n + K) instead of O((n+W) log(n+W)) — no sorting of
    the 193k WT array per call. NC counts are small non-negative ints,
    so this is exact (verified np.allclose vs scipy)."""
    a = np.asarray(a)
    b = np.asarray(b)
    K = int(max(a.max(), b.max()))
    ca = np.cumsum(np.bincount(a, minlength=K + 1)[:K + 1]) / a.size
    cb = np.cumsum(np.bincount(b, minlength=K + 1)[:K + 1]) / b.size
    return float(np.sum(np.abs(ca - cb)))


def nc_shift_null(wt_nc_samples, n, n_iter, rng):
    """Null distribution of the variant-vs-WT NC Wasserstein-1 statistic
    under H0 ('variant behaves like WT'): draw n WT NC values WITH
    replacement (NOT NC-matched — the point is to detect an NC shift)
    and measure Wasserstein-1 vs the full WT NC distribution. Depends
    only on n, so it is computed once per stratum and reused.

    The WT CDF is fixed across all iterations; precompute it once and
    score each resample with the exact integer W1 (was 10000x
    scipy.wasserstein_distance over the 193k WT array per stratum — the
    dominant production cost). rng.choice draws are unchanged, so the
    null values are identical to the scipy version.
    """
    wt = np.asarray(wt_nc_samples)
    K = int(wt.max())
    cb = np.cumsum(np.bincount(wt, minlength=K + 1)[:K + 1]) / wt.size
    null_w = np.empty(n_iter, dtype=np.float64)
    for i in range(n_iter):
        s = rng.choice(wt_nc_samples, size=n, replace=True)
        ca = np.cumsum(np.bincount(s, minlength=K + 1)[:K + 1]) / s.size
        null_w[i] = np.sum(np.abs(ca - cb))
    return null_w


def nc_shift_stats(var_nc_samples, wt_nc_samples, null_w):
    """Per-variant nucleosome-count shift readout (distinct from the
    per-position nuc-track occupancy test). Signed ΔNC mean for
    direction/magnitude + Wasserstein-1 with an empirical p vs the
    shared null. obs uses the same exact integer W1 as the null so the
    null_w >= obs_w comparison is on an identical metric implementation.
    """
    if len(var_nc_samples) == 0 or len(wt_nc_samples) == 0:
        return None
    obs_w = _w1_nonneg_int(var_nc_samples, wt_nc_samples)
    mv = float(np.mean(var_nc_samples))
    mw = float(np.mean(wt_nc_samples))
    p = float((np.sum(null_w >= obs_w) + 1) / (len(null_w) + 1))
    return {'nc_mean_variant': mv, 'nc_mean_wt': mw,
            'nc_delta': mv - mw, 'nc_wasserstein': obs_w,
            'nc_shift_p': p}


def build_strata(rd, variant_groups, min_nuc=None, max_nuc=None,
                 n_tol=0.10, nc_dist=0.30, stratify=True):
    """Group testable variants whose per-variant null would be ~identical
    (similar N and NC distribution) so one null can be reused.

    Greedy, sorted by N. A variant joins a stratum if its N is within
    n_tol (relative) of the stratum AND its NC distribution is within
    Wasserstein-1 nc_dist of the stratum representative (fixed at
    creation, like the centroid-fixed barcode clustering). The stratum
    null depth is the MIN member N (conservative: the null is at least
    as wide as any member's true sampling noise). stratify=False puts
    every variant in its own stratum (== independent per-variant null).

    Returns list of {rep_n, rep_nc, members:[(vid, var_f, var_n, vnc)]}.
    """
    entries = []
    for vid, var_indices in variant_groups.items():
        var_f = rd.get_indices_filtered(var_indices, min_nuc=min_nuc,
                                        max_nuc=max_nuc)
        var_n = len(var_f)
        if var_n < 20:
            continue
        vnc = compute_ground_truth_nc(rd, var_f, min_nuc=min_nuc,
                                      max_nuc=max_nuc)
        if vnc is None:
            continue
        entries.append((vid, var_f, var_n, vnc))
    entries.sort(key=lambda e: e[2])

    strata = []
    for vid, var_f, var_n, vnc in entries:
        placed = False
        if stratify:
            for st in strata:
                if (abs(var_n - st['rep_n']) <= n_tol * max(1, st['rep_n'])
                        and _nc_wasserstein(vnc, st['rep_nc']) <= nc_dist):
                    st['members'].append((vid, var_f, var_n, vnc))
                    if var_n < st['rep_n']:
                        st['rep_n'] = var_n  # conservative: widest null
                    placed = True
                    break
        if not placed:
            strata.append({'rep_n': var_n, 'rep_nc': vnc,
                           'members': [(vid, var_f, var_n, vnc)]})
    return strata


def run_variant_testing_parallel(rd, wt_idx, variant_groups,
                                  analysis_region,
                                  promoter_start, promoter_end,
                                  min_nuc=None, max_nuc=None,
                                  random_seed=42,
                                  n_null_iterations=2000,
                                  cluster_threshold_quantile=0.95,
                                  absolute_delta_threshold=None,
                                  gap_tolerance=2, merge_distance=5,
                                  n_workers=1, stratify=True,
                                  n_tol=0.10, nc_dist=0.30,
                                  mde_alpha=0.05):
    """Test all variants. Nulls are built once per stratum (variants
    with ~identical N and NC distribution share a null — the compute
    lever) and reused across members. Each variant's observed delta
    uses its OWN full read set (Option-A reference from the stratum
    null). The expensive null build is parallelized across n_workers
    by run_null_calibration's shared-memory path.
    """
    a_start, a_end = analysis_region
    analysis_length = a_end - a_start
    prom_rel_s = max(0, promoter_start - a_start)
    prom_rel_e = min(analysis_length, promoter_end - a_start)
    prom_slice = slice(prom_rel_s, prom_rel_e)
    bin_labels = rd.bin_labels

    strata = build_strata(rd, variant_groups, min_nuc=min_nuc,
                          max_nuc=max_nuc, n_tol=n_tol, nc_dist=nc_dist,
                          stratify=stratify)
    n_members = sum(len(s['members']) for s in strata)
    logging.info(f"  Null stratification: {n_members} testable variants "
                 f"-> {len(strata)} strata (stratify={stratify}, "
                 f"n_tol={n_tol}, nc_dist={nc_dist})")

    needs_col_slice = (rd.analysis_length != analysis_length
                       or rd.analysis_start != a_start)
    col_indices = (np.arange(a_start, a_end) if needs_col_slice else None)

    # WT NC samples (full WT, not NC-matched) for the NC-shift readout.
    wt_f_nc = rd.get_indices_filtered(wt_idx, min_nuc=min_nuc,
                                      max_nuc=max_nuc)
    wt_nc_samples = rd.nuc_counts[wt_f_nc]
    wt_nc_samples = wt_nc_samples[wt_nc_samples >= 0]

    # PERF: the WT slice / NC pools / shared memory and the worker pool
    # are INVARIANT across variants — build them ONCE here and reuse for
    # every stratum, instead of run_null_calibration rebuilding them
    # 149x. Numerically identical (same arrays, same seeds, same
    # chunking); this is purely "build once, reuse" vs "rebuild each".
    use_pool = n_workers > 1 and n_null_iterations > 1
    wt_ctx = build_wt_null_context(rd, wt_idx, min_nuc, max_nuc,
                                   analysis_region, make_shared=use_pool)
    persistent_ex = (ProcessPoolExecutor(max_workers=n_workers)
                     if use_pool else None)
    if persistent_ex is not None:
        # Pre-spawn all workers now so the FIRST stratum does not pay
        # the one-time cold-pool spawn (observed ~minutes at 64 workers
        # on the cluster vs ~seconds warm). A tiny sleep per task makes
        # each task land on a distinct worker so all n_workers spin up.
        # Numerically inert — no statistics touched.
        t_warm = time.time()
        list(persistent_ex.map(_pool_warmup, range(n_workers)))
        logging.info(f"  Worker pool pre-warmed ({n_workers} workers, "
                     f"{time.time() - t_warm:.1f}s)")

    all_results = []
    try:
      for si, st in enumerate(tqdm(strata, desc="Strata",
                                   disable=not HAS_TQDM)):
        null_res = run_null_calibration(
            rd, wt_idx, st['rep_nc'], coverage_level=st['rep_n'],
            min_nuc=min_nuc, max_nuc=max_nuc,
            n_iterations=n_null_iterations,
            random_seed=stable_variant_seed(random_seed, f"stratum_{si}"),
            analysis_region=analysis_region,
            cluster_threshold_quantile=cluster_threshold_quantile,
            absolute_delta_threshold=absolute_delta_threshold,
            gap_tolerance=gap_tolerance, merge_distance=merge_distance,
            n_workers=n_workers, wt_ctx=wt_ctx, executor=persistent_ex)
        # NC-shift null depends only on N -> one per stratum (rep_n).
        nc_null = (nc_shift_null(
            wt_nc_samples, st['rep_n'], n_null_iterations,
            np.random.default_rng(
                stable_variant_seed(random_seed, f"ncshift_{si}")))
            if len(wt_nc_samples) else np.zeros(n_null_iterations))
        for vid, var_f, var_n, vnc in st['members']:
            if needs_col_slice:
                def var_occ_fn(label, _vf=var_f):
                    return rd.coverage_matrices[label][
                        np.ix_(_vf, col_indices)].mean(axis=0)
            else:
                def var_occ_fn(label, _vf=var_f):
                    return rd.coverage_matrices[label][_vf].mean(axis=0)
            vr = _compute_variant_result(
                vid, variant_groups[vid], var_f, var_n, var_occ_fn,
                null_res, var_n, bin_labels, a_start, analysis_length,
                prom_rel_s, prom_rel_e, prom_slice,
                cluster_threshold_quantile, absolute_delta_threshold,
                gap_tolerance, merge_distance, mde_alpha=mde_alpha)
            if vr is not None:
                vnc_samp = rd.nuc_counts[var_f]
                vnc_samp = vnc_samp[vnc_samp >= 0]
                ncs = nc_shift_stats(vnc_samp, wt_nc_samples, nc_null)
                if ncs is not None:
                    vr.update(ncs)
                all_results.append(vr)
    finally:
        if persistent_ex is not None:
            persistent_ex.shutdown()
        if wt_ctx.get('shm_objects'):
            cleanup_shared_memory(wt_ctx['shm_objects'])
    logging.info(f"  Tested {len(all_results)} variants "
                 f"({len(strata)} null builds)")
    return all_results


# ============================================================================
# Phase 2 — TF binding-motif layer (plan §2.3)
#
# Implemented as a post-hoc annotation + cross-variant aggregation layer on
# top of the per-variant results. It deliberately does NOT touch the
# validated Phase 1 statistical hot path (run_null_calibration /
# _compute_variant_result); it only reads delta_obs and the already-detected
# significant clusters and the calibrated per-variant FDR.
# ============================================================================

DEFAULT_MOTIF_CFG = {
    'bins': ('TF', 'sub_TF'),
    'direction': 'loss',          # 'loss' | 'gain' | 'both'
    'sign_consistency': 0.90,
    'min_width': 5,
    'max_width': 25,
    'require_variant_overlap': False,  # per-variant: annotate only
    'fdr': 0.10,                  # variant gate: variant_fdr_q < fdr
    'density_threshold': 2,       # >= N distinct variants -> motif call
}


def build_motif_cfg(bins=None, direction='loss', sign_consistency=0.90,
                     min_width=5, max_width=25,
                     require_variant_overlap=False, fdr=0.10,
                     density_threshold=2):
    """Assemble the motif-layer config from CLI/Snakefile values, with the
    plan §2.3 defaults. Kept as a plain dict so it pickles trivially."""
    return {
        'bins': tuple(bins) if bins else DEFAULT_MOTIF_CFG['bins'],
        'direction': direction,
        'sign_consistency': float(sign_consistency),
        'min_width': int(min_width),
        'max_width': int(max_width),
        'require_variant_overlap': bool(require_variant_overlap),
        'fdr': float(fdr),
        'density_threshold': int(density_threshold),
    }


def load_reference_sequence(fasta_path, ref_name):
    """Load one contig's sequence (uppercase str) from a FASTA, or None.

    Used only for motif DNA extraction (§4 decision #1: Stage 3 takes an
    optional --reference rather than Stage 2 carrying the sequence)."""
    if not fasta_path:
        return None
    try:
        fa = pysam.FastaFile(fasta_path)
    except Exception as e:
        logging.warning(f"Could not open --reference {fasta_path}: {e}")
        return None
    try:
        if ref_name not in fa.references:
            logging.warning(
                f"--reference has no contig '{ref_name}' "
                f"(has: {list(fa.references)[:5]}...); "
                f"motif DNA extraction disabled")
            return None
        seq = fa.fetch(ref_name).upper()
        logging.info(f"Reference loaded for motifs: {ref_name} "
                     f"({len(seq):,} bp)")
        return seq
    finally:
        fa.close()


def _cluster_sign_consistency(delta_obs, c):
    """Fraction of positions inside the cluster whose signed Δ matches the
    cluster's dominant sign. delta_obs is the analysis-window array; the
    cluster's start/end index into it (relative coords)."""
    seg = np.asarray(delta_obs[c['start']:c['end'] + 1], dtype=np.float64)
    if seg.size == 0:
        return 0.0
    dom = 1.0 if c.get('mean_signed_delta', 0.0) >= 0 else -1.0
    return float(np.mean(np.sign(seg) == dom))


def annotate_variant_motifs(all_variant_results, ref_seq, analysis_region,
                            motif_cfg):
    """Annotate every per-variant cluster in place with the motif-layer
    fields and flag motif clusters per plan §2.3.

    Adds to each cluster dict: sign_consistency, ref_sequence (if
    ref_seq), variant_distance (causal-variant gap; 0 = overlap),
    variant_overlap (bool), is_motif (bool).

    Coordinate note: variant IDs are 1-based ref positions
    (02_call_variants.py emits position = 0-based + 1); cluster
    abs_start/abs_end are 0-based ref positions. We convert the variant
    to 0-based before comparing/extracting.
    """
    a_start, a_end = analysis_region
    bins = set(motif_cfg['bins'])
    want_dir = motif_cfg['direction']
    sc_min = motif_cfg['sign_consistency']
    w_min, w_max = motif_cfg['min_width'], motif_cfg['max_width']
    req_ov = motif_cfg['require_variant_overlap']
    n_ref = len(ref_seq) if ref_seq else 0
    ref_base_mismatches = 0
    ref_base_checked = 0

    for vr in all_variant_results:
        vid = vr['variant_id']
        pos1, vref, valt, vct = parse_variant_id_fields(vid)
        var_pos0 = (pos1 - 1) if pos1 is not None else None
        vr['variant_pos0'] = var_pos0
        vr['variant_alt'] = valt
        vr['variant_change_type'] = vct
        # Empirical coordinate self-check: ref base at the SNV should
        # match the variant's recorded ref base (locks 1-based vs
        # 0-based convention; logged, not fatal).
        if (ref_seq and var_pos0 is not None and vct == 'snv'
                and vref and 0 <= var_pos0 < n_ref):
            ref_base_checked += 1
            if ref_seq[var_pos0].upper() != vref.upper():
                ref_base_mismatches += 1

        for label in list(vr.keys()):
            lr = vr.get(label)
            if not isinstance(lr, dict) or 'delta_obs' not in lr:
                continue
            delta_obs = lr['delta_obs']
            seen = set()
            for clist_key in ('significant_clusters', 'all_promoter_clusters'):
                for c in lr.get(clist_key, []):
                    cid = id(c)
                    if cid in seen:
                        continue
                    seen.add(cid)
                    c['sign_consistency'] = _cluster_sign_consistency(
                        delta_obs, c)
                    a_s = c.get('abs_start', c['start'] + a_start)
                    a_e = c.get('abs_end', c['end'] + a_start)
                    c['abs_start'] = a_s
                    c['abs_end'] = a_e
                    if ref_seq and 0 <= a_s <= a_e < n_ref:
                        c['ref_sequence'] = ref_seq[a_s:a_e + 1]
                    else:
                        c['ref_sequence'] = ''
                    if var_pos0 is None:
                        c['variant_overlap'] = False
                        c['variant_distance'] = -999999
                    elif a_s <= var_pos0 <= a_e:
                        c['variant_overlap'] = True
                        c['variant_distance'] = 0
                    else:
                        c['variant_overlap'] = False
                        c['variant_distance'] = int(
                            var_pos0 - a_s if var_pos0 < a_s
                            else var_pos0 - a_e)
                    dir_ok = (want_dir == 'both'
                              or c.get('direction') == want_dir)
                    c['is_motif'] = bool(
                        label in bins
                        and w_min <= c['width'] <= w_max
                        and c['sign_consistency'] >= sc_min
                        and dir_ok
                        and (not req_ov or c['variant_overlap']))

    if ref_base_checked:
        frac = ref_base_mismatches / ref_base_checked
        msg = (f"Motif coord check: {ref_base_mismatches}/"
               f"{ref_base_checked} SNV ref bases disagree with "
               f"--reference ({frac:.1%})")
        if frac > 0.05:
            logging.warning(
                msg + " — possible coordinate/contig mismatch; "
                "motif DNA may be off by one or wrong contig")
        else:
            logging.info(msg)


def aggregate_motifs(all_variant_results, ref_seq, ref_name,
                     analysis_region, motif_cfg):
    """Cross-variant aggregation — the primary biological deliverable
    (plan §2.3). For each motif bin and direction, build a
    disruption-density track (# distinct FDR-significant variants whose
    sign-consistent motif cluster covers each reference position), call
    motif intervals where density >= threshold, and for each motif emit
    the reference DNA plus a per-position/per-base sensitivity profile.

    Returns {'tracks': {(bin,dir): np.int32[analysis_length]},
             'motifs': [ {...} ], 'analysis_start', 'analysis_length',
             'ref_name'} .
    """
    a_start, a_end = analysis_region
    analysis_length = a_end - a_start
    bins = list(motif_cfg['bins'])
    if motif_cfg['direction'] == 'both':
        directions = ['loss', 'gain']
    else:
        directions = [motif_cfg['direction']]
    fdr = motif_cfg['fdr']
    dens_thr = motif_cfg['density_threshold']
    w_min = motif_cfg['min_width']
    req_ov = motif_cfg['require_variant_overlap']
    n_ref = len(ref_seq) if ref_seq else 0
    BASES = ('A', 'C', 'G', 'T')
    bidx = {b: i for i, b in enumerate(BASES)}

    sig_variants = [vr for vr in all_variant_results
                    if vr.get('variant_fdr_q', 1.0) < fdr]

    tracks = {}
    motifs = []
    for label in bins:
        for d in directions:
            density = np.zeros(analysis_length, dtype=np.int32)
            # per-variant: union mask of its qualifying motif clusters,
            # plus the clusters kept for the overlap/sensitivity step.
            per_var = []
            for vr in sig_variants:
                lr = vr.get(label)
                if not isinstance(lr, dict):
                    continue
                qcl = [c for c in lr.get('significant_clusters', [])
                       if c.get('is_motif')
                       and (motif_cfg['direction'] == 'both'
                            or c.get('direction') == d)
                       and (not req_ov or c.get('variant_overlap'))]
                if not qcl:
                    continue
                mask = np.zeros(analysis_length, dtype=bool)
                for c in qcl:
                    mask[c['start']:c['end'] + 1] = True
                density += mask.astype(np.int32)
                per_var.append((vr, qcl, mask))
            tracks[(label, d)] = density

            # Call motif intervals: contiguous runs >= threshold.
            hot = density >= dens_thr
            i = 0
            while i < analysis_length:
                if not hot[i]:
                    i += 1
                    continue
                j = i
                while j < analysis_length and hot[j]:
                    j += 1
                m_s, m_e = i, j - 1          # relative, inclusive
                i = j
                if (m_e - m_s + 1) < w_min:
                    continue
                abs_s, abs_e = m_s + a_start, m_e + a_start
                L = m_e - m_s + 1
                if ref_seq and 0 <= abs_s <= abs_e < n_ref:
                    mseq = ref_seq[abs_s:abs_e + 1]
                else:
                    mseq = ''
                sens_count = np.zeros((L, 4), dtype=np.int32)
                sens_effect = np.zeros((L, 4), dtype=np.float64)
                contributors = []
                for vr, qcl, mask in per_var:
                    if not mask[m_s:m_e + 1].any():
                        continue
                    contributors.append(vr['variant_id'])
                    vp0 = vr.get('variant_pos0')
                    valt = (vr.get('variant_alt') or '').upper()
                    if (vr.get('variant_change_type') == 'snv'
                            and vp0 is not None
                            and abs_s <= vp0 <= abs_e
                            and valt in bidx):
                        # strongest overlapping motif cluster's signed Δ
                        best = None
                        for c in qcl:
                            if c['end'] >= m_s and c['start'] <= m_e:
                                msd = c.get('mean_signed_delta', 0.0)
                                if best is None or abs(msd) > abs(best):
                                    best = msd
                        r = vp0 - abs_s
                        b = bidx[valt]
                        sens_count[r, b] += 1
                        if best is not None:
                            sens_effect[r, b] += best
                with np.errstate(invalid='ignore', divide='ignore'):
                    sens_mean = np.where(
                        sens_count > 0,
                        sens_effect / np.maximum(sens_count, 1), 0.0)
                motifs.append({
                    'bin': label, 'direction': d,
                    'abs_start': int(abs_s), 'abs_end': int(abs_e),
                    'width': int(L), 'ref_sequence': mseq,
                    'n_variants': len(set(contributors)),
                    'peak_density': int(density[m_s:m_e + 1].max()),
                    'contributing_variant_ids': sorted(set(contributors)),
                    'sensitivity_count': sens_count,
                    'sensitivity_mean_signed_delta': sens_mean,
                    'base_order': ''.join(BASES),
                })

    motifs.sort(key=lambda m: (-m['n_variants'], -m['peak_density']))
    logging.info(
        f"  Motif aggregation: {len(sig_variants)} FDR<{fdr} variants "
        f"-> {len(motifs)} motif call(s) "
        f"(bins={bins}, dir={motif_cfg['direction']}, "
        f"density>={dens_thr})")
    return {'tracks': tracks, 'motifs': motifs,
            'analysis_start': a_start,
            'analysis_length': analysis_length,
            'ref_name': ref_name}


# ============================================================================
# Phase 3 — protein co-occupancy dependency (plan §10 / §2.4–2.5)
#
# P3.2: site-pair construction. Operates ONLY on the in-memory Phase-2
# results (all_variant_results, motif_result) + rd. Does not touch any
# validated Phase-1/2 path. Default OFF (--enable-co-occupancy).
# ============================================================================

DEFAULT_COOCCUPANCY_CFG = {
    'site_occupancy_frac': 0.5,    # read "occupies" a site if >= this
                                   # fraction of site positions footprinted
    'site_merge_overlap': 0.5,     # merge same-bin candidate intervals with
                                   # reciprocal overlap >= this -> 1 canonical
    'min_site_separation': 30,     # bp gap required between site1 and site2
    'min_variant_instruments': 3,  # site1 needs >= this many independent
                                   # disrupting SNVs to be testable
    'motif_fdr': 0.10,             # variant gate (same as the motif layer)
    'min_stratum': 25,             # min WT reads in each O1 stratum
    'call_q': 0.10,                # BH q threshold for a co-occ CALL
    'min_consistency': 0.70,       # frac of instruments agreeing in sign
                                   # required to CALL a dependency
                                   # (Phase-2 sign-consistency analogue:
                                   # cooperativity = independent SNVs
                                   # AGREE on direction, not just a
                                   # nonzero mean excess)
    'secondary_screen': True,      # §2.5 control: validate calls by
                                   # re-testing on reads with NO
                                   # competing site-2-local raw SNV
    'secondary_window': 10,        # bp pad around site 2 for "local"
                                   # secondary raw mutations
}


def build_cooccupancy_cfg(site_occupancy_frac=0.5, site_merge_overlap=0.5,
                          min_site_separation=30,
                          min_variant_instruments=3, motif_fdr=0.10,
                          min_stratum=25, call_q=0.10,
                          min_consistency=0.70, secondary_screen=True,
                          secondary_window=10):
    return {
        'site_occupancy_frac': float(site_occupancy_frac),
        'site_merge_overlap': float(site_merge_overlap),
        'min_site_separation': int(min_site_separation),
        'min_variant_instruments': int(min_variant_instruments),
        'motif_fdr': float(motif_fdr),
        'min_stratum': int(min_stratum),
        'call_q': float(call_q),
        'min_consistency': float(min_consistency),
        'secondary_screen': bool(secondary_screen),
        'secondary_window': int(secondary_window),
    }


def _recip_overlap(a_s, a_e, b_s, b_e):
    """Reciprocal overlap fraction of two inclusive intervals:
    overlap_len / min(len_a, len_b). 0 if disjoint."""
    ov = min(a_e, b_e) - max(a_s, b_s) + 1
    if ov <= 0:
        return 0.0
    return ov / min(a_e - a_s + 1, b_e - b_s + 1)


def _merge_bin_intervals(intervals, merge_overlap):
    """Greedily collapse a bin's candidate intervals: sort by start, fuse
    into a running interval while reciprocal overlap >= merge_overlap.
    Returns list of (start, end) canonical intervals (union extents)."""
    if not intervals:
        return []
    iv = sorted(intervals)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        cs, ce = out[-1]
        if _recip_overlap(cs, ce, s, e) >= merge_overlap:
            out[-1][1] = max(ce, e)
            out[-1][0] = min(cs, s)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def build_canonical_sites(all_variant_results, motif_result, bin_labels,
                          cfg):
    """Canonical footprint sites = (bin, abs_start, abs_end). Candidate
    intervals are pooled from the cross-variant motif calls AND every
    FDR-significant variant's significant clusters (the decided site-2
    set), then collapsed per bin so the many per-variant clusters over
    one element become a single site. Returns a list of dicts with a
    stable site_id and whether the site is also a called motif.
    """
    fdr = cfg['motif_fdr']
    by_bin = {b: [] for b in bin_labels}
    motif_intervals = set()
    if motif_result is not None:
        for m in motif_result['motifs']:
            by_bin.setdefault(m['bin'], []).append(
                (m['abs_start'], m['abs_end']))
            motif_intervals.add((m['bin'], m['abs_start'], m['abs_end']))
    for vr in all_variant_results:
        if vr.get('variant_fdr_q', 1.0) >= fdr:
            continue
        for label in bin_labels:
            lr = vr.get(label)
            if not isinstance(lr, dict):
                continue
            for c in lr.get('significant_clusters', []):
                by_bin.setdefault(label, []).append(
                    (c['abs_start'], c['abs_end']))

    sites = []
    for label in bin_labels:
        for (s, e) in _merge_bin_intervals(by_bin.get(label, []),
                                           cfg['site_merge_overlap']):
            is_motif = any(
                mb == label and _recip_overlap(s, e, ms, me) > 0
                for (mb, ms, me) in motif_intervals)
            sites.append({
                'site_id': f"{label}:{s}-{e}",
                'bin': label, 'abs_start': int(s), 'abs_end': int(e),
                'is_motif_site': bool(is_motif),
            })
    return sites


def _site_columns(site, analysis_start, analysis_length):
    """Coverage-matrix column indices for a site (abs ref coords ->
    analysis-window columns), clipped to the window."""
    c0 = max(0, site['abs_start'] - analysis_start)
    c1 = min(analysis_length - 1, site['abs_end'] - analysis_start)
    if c1 < c0:
        return np.empty(0, dtype=np.int64)
    return np.arange(c0, c1 + 1, dtype=np.int64)


def read_site_occupied(rd, site, read_idx, occ_frac):
    """Per-read binary state at a site: True if the read's footprint in
    the site's bin covers >= occ_frac of the site's positions. Returns a
    bool array aligned to read_idx."""
    cols = _site_columns(site, rd.analysis_start, rd.analysis_length)
    read_idx = np.asarray(read_idx)
    if cols.size == 0 or read_idx.size == 0:
        return np.zeros(read_idx.size, dtype=bool)
    sub = rd.coverage_matrices[site['bin']][np.ix_(read_idx, cols)]
    return sub.mean(axis=1) >= occ_frac


def assign_variant_site1(vr, sites, bin_labels, cfg):
    """The site this variant disrupts: among canonical sites that
    contain the causal SNV and overlap one of the variant's is_motif
    significant clusters, pick the one with the largest such cluster
    (Σ|Δ|). Returns the site dict or None (variant not a usable
    instrument)."""
    if vr.get('variant_fdr_q', 1.0) >= cfg['motif_fdr']:
        return None
    vp = vr.get('variant_pos0')
    if vp is None:
        return None
    best = None
    best_sum = -1.0
    for label in bin_labels:
        lr = vr.get(label)
        if not isinstance(lr, dict):
            continue
        for c in lr.get('significant_clusters', []):
            if not c.get('is_motif'):
                continue
            for site in sites:
                if site['bin'] != label:
                    continue
                if not (site['abs_start'] <= vp <= site['abs_end']):
                    continue
                if _recip_overlap(site['abs_start'], site['abs_end'],
                                  c['abs_start'], c['abs_end']) <= 0:
                    continue
                if c['sum_abs_delta'] > best_sum:
                    best_sum = c['sum_abs_delta']
                    best = site
    return best


def build_site_pairs(all_variant_results, motif_result, rd, bin_labels,
                     cfg):
    """Assemble co-occupancy testing units (P3.2 output):

      sites              : canonical sites (build_canonical_sites)
      site1_instruments  : {site_id: [variant_id, ...]} — independent
                           SNVs that disrupt that site (the "instruments")
      pairs              : [(site1_id, site2_id), ...] for every site1
                           with >= min_variant_instruments, paired with
                           each canonical site2 that is >= min_site_
                           separation bp away (the decided site-2 set)

    No statistics here — just the structure P3.3 will test.
    """
    sites = build_canonical_sites(all_variant_results, motif_result,
                                  bin_labels, cfg)
    by_id = {s['site_id']: s for s in sites}
    instruments = defaultdict(list)
    for vr in all_variant_results:
        s1 = assign_variant_site1(vr, sites, bin_labels, cfg)
        if s1 is not None:
            instruments[s1['site_id']].append(vr['variant_id'])

    sep = cfg['min_site_separation']
    pairs = []
    for s1_id, vids in instruments.items():
        if len(vids) < cfg['min_variant_instruments']:
            continue
        s1 = by_id[s1_id]
        for s2 in sites:
            if s2['site_id'] == s1_id:
                continue
            gap = max(s1['abs_start'] - s2['abs_end'],
                      s2['abs_start'] - s1['abs_end'])
            if gap < sep:           # overlapping/too close -> not distal
                continue
            pairs.append((s1_id, s2['site_id']))

    logging.info(
        f"  Co-occupancy P3.2: {len(sites)} canonical sites, "
        f"{sum(1 for v in instruments.values() if len(v) >= cfg['min_variant_instruments'])}"
        f" testable site1 (>= {cfg['min_variant_instruments']} instruments), "
        f"{len(pairs)} (site1,site2) pairs")
    return {'sites': sites, 'sites_by_id': by_id,
            'site1_instruments': dict(instruments),
            'pairs': pairs, 'cfg': cfg}


# ---- P3.3: WT conditional model + per-variant excess test --------------
#
# Statistics are factored into pure array-level helpers (unit-testable on
# synthetic ground truth) with thin rd-facing wrappers. Null = conditional
# resampling with the WT reads BOOTSTRAPPED each iteration so it carries
# the WT-conditional's own estimation uncertainty (Phase-1 lesson: plug-in
# nulls are anti-conservative). The WT-conditional bootstrap is computed
# ONCE per site-pair and reused across all that site-1's instruments.

def _wt_conditional_from_states(o1w, o2w, n_iter, rng, min_stratum):
    """WT P(O2|O1) point estimates + a bootstrap distribution.

    The WT joint is fully summarized by the 4 cell counts
    (O1,O2)∈{0,1}², so a WT read-resample is exactly a multinomial
    draw over those 4 cells — O(n_iter), not O(n_iter·W), and
    identical in distribution to resampling reads. Returns None if
    either O1 stratum has < min_stratum WT reads (pair untestable).
    """
    o1w = np.asarray(o1w, dtype=bool)
    o2w = np.asarray(o2w, dtype=bool)
    N = o1w.size
    c11 = int(np.count_nonzero(o1w & o2w))
    c10 = int(np.count_nonzero(o1w & ~o2w))
    c01 = int(np.count_nonzero(~o1w & o2w))
    c00 = int(np.count_nonzero(~o1w & ~o2w))
    n1 = c10 + c11          # WT reads with O1=1
    n0 = c00 + c01          # WT reads with O1=0
    if n1 < min_stratum or n0 < min_stratum or N == 0:
        return None
    p1 = c11 / n1
    p0 = c01 / n0
    probs = np.array([c00, c01, c10, c11], dtype=np.float64) / N
    m = rng.multinomial(N, probs, size=n_iter)        # (B,4): 00,01,10,11
    m00, m01, m10, m11 = m[:, 0], m[:, 1], m[:, 2], m[:, 3]
    d1 = m10 + m11
    d0 = m00 + m01
    overall = (m01 + m11) / N                          # fallback if a
    p1_boot = np.where(d1 > 0, m11 / np.maximum(d1, 1), overall)
    p0_boot = np.where(d0 > 0, m01 / np.maximum(d0, 1), overall)
    return {'p1': p1, 'p0': p0, 'p1_boot': p1_boot,
            'p0_boot': p0_boot, 'n_w1': n1, 'n_w0': n0,
            'counts': (c00, c01, c10, c11)}


def _cooccupancy_test_from_states(o1v, o2v, wt_cond, rng,
                                  chunk=2_000_000):
    """Per-variant excess test given the variant's per-read site states
    and a precomputed WT-conditional bootstrap (reused across the
    site-1's instruments). H0: site 2 follows WT's P(O2|O1) given the
    variant's realized O1. Two-sided empirical p; signed excess vs the
    parametric prediction as the effect-size readout."""
    o1v = np.asarray(o1v, dtype=bool)
    o2v = np.asarray(o2v, dtype=bool)
    n_v = o1v.size
    if n_v == 0:
        return None
    p1, p0 = wt_cond['p1'], wt_cond['p0']
    p1b, p0b = wt_cond['p1_boot'], wt_cond['p0_boot']
    B = p1b.size
    obs = float(o2v.mean())
    # Parametric prediction (effect-size readout): expected site-2
    # occupancy if site 2 only responds via the WT O2|O1 channel.
    predicted = float(np.where(o1v, p1, p0).mean())

    # Null: per bootstrap b, draw each variant read's O2 ~
    # Bernoulli(p_{O1v,i}^*(b)); statistic = mean(O2). Chunk over B so
    # the (B × n_v) draw stays bounded in RAM.
    null_stats = np.empty(B, dtype=np.float64)
    step = max(1, chunk // max(n_v, 1))
    for s in range(0, B, step):
        e = min(B, s + step)
        P = np.where(o1v[None, :],
                     p1b[s:e, None], p0b[s:e, None])      # (k, n_v)
        U = rng.random((e - s, n_v))
        null_stats[s:e] = (U < P).mean(axis=1)

    m = float(null_stats.mean())
    p_two = float((1 + np.count_nonzero(
        np.abs(null_stats - m) >= abs(obs - m))) / (B + 1))
    excess = obs - predicted
    return {
        'n_var': int(n_v),
        'obs_site2_occ': obs,
        'predicted_site2_occ': predicted,
        'excess': float(excess),
        'p_two_sided': p_two,
        'direction': 'site2_loss' if excess < 0 else 'site2_gain',
        'null_mean': m,
        'null_std': float(null_stats.std()),
    }


def wt_conditional_bootstrap(rd, site1, site2, wt_idx, occ_frac,
                             n_iter, rng, min_stratum):
    """rd-facing: extract WT per-read site states then build the
    reusable WT-conditional bootstrap. None -> pair untestable."""
    o1w = read_site_occupied(rd, site1, wt_idx, occ_frac)
    o2w = read_site_occupied(rd, site2, wt_idx, occ_frac)
    return _wt_conditional_from_states(o1w, o2w, n_iter, rng,
                                       min_stratum)


def cooccupancy_variant_test(rd, site1, site2, var_idx, occ_frac,
                             wt_cond, rng):
    """rd-facing per-variant test (uses the pair's precomputed
    wt_cond)."""
    o1v = read_site_occupied(rd, site1, var_idx, occ_frac)
    o2v = read_site_occupied(rd, site2, var_idx, occ_frac)
    return _cooccupancy_test_from_states(o1v, o2v, wt_cond, rng)


# ---- P3.4: cross-variant aggregation (the PRIMARY evidence) + BH ------
#
# A real cooperative dependency is corroborated by MANY INDEPENDENT SNVs
# that all disrupt site 1 showing a consistent, directional site-2 shift
# beyond the WT channel. We test a pair-level weighted-mean-excess
# statistic against a JOINT null: per bootstrap b, every instrument's H0
# site-2 mean is drawn under the SAME WT-conditional bootstrap b (shared
# WT-estimation noise -> the instruments' nulls are correlated; the joint
# resample captures that correlation, which a per-variant-then-combine
# scheme would miss — the Phase-1 calibration lesson). BH across pairs.

def _cooccupancy_aggregate_from_states(inst_states, wt_cond, rng,
                                       weight='equal',
                                       chunk=2_000_000):
    """inst_states: list of (o1_k, o2_k) per instrument (independent
    SNVs disrupting the shared site 1). Returns the pair-level
    aggregate test + per-instrument readouts. Weight 'equal' (each
    independent instrument = one evidence unit; robust, the truer
    'many SNVs agree' test) or 'n' (read-count weighted)."""
    p1, p0 = wt_cond['p1'], wt_cond['p0']
    p1b, p0b = wt_cond['p1_boot'], wt_cond['p0_boot']
    B = p1b.size
    K = len(inst_states)
    if K == 0:
        return None
    obs_k = np.empty(K)
    pred_k = np.empty(K)
    w_k = np.empty(K)
    per_inst = []
    null_T = np.zeros(B, dtype=np.float64)
    for k, (o1, o2) in enumerate(inst_states):
        o1 = np.asarray(o1, dtype=bool)
        o2 = np.asarray(o2, dtype=bool)
        n_k = o1.size
        obs_k[k] = o2.mean() if n_k else 0.0
        pred_k[k] = np.where(o1, p1, p0).mean() if n_k else 0.0
        w_k[k] = n_k if weight == 'n' else 1.0
        # this instrument's H0 null mean(O2) for every bootstrap b
        nb = np.empty(B)
        step = max(1, chunk // max(n_k, 1))
        for s in range(0, B, step):
            e = min(B, s + step)
            P = np.where(o1[None, :], p1b[s:e, None], p0b[s:e, None])
            U = rng.random((e - s, n_k))
            nb[s:e] = (U < P).mean(axis=1)
        # accumulate into the JOINT pair statistic (same b across k)
        null_T += w_k[k] * (nb - pred_k[k])
        per_inst.append({'n': int(n_k), 'obs': float(obs_k[k]),
                         'predicted': float(pred_k[k]),
                         'excess': float(obs_k[k] - pred_k[k])})
    W = w_k.sum()
    excess_k = obs_k - pred_k
    T_obs = float((w_k * excess_k).sum() / W)
    null_T /= W
    m = float(null_T.mean())
    p_two = float((1 + np.count_nonzero(
        np.abs(null_T - m) >= abs(T_obs - m))) / (B + 1))
    # one-sided for the primary hypothesis (site-2 LOSS beyond channel)
    p_loss = float((1 + np.count_nonzero(null_T <= T_obs)) / (B + 1))
    n_consistent = int(np.count_nonzero(
        np.sign(excess_k) == np.sign(T_obs))) if T_obs != 0 else 0
    return {
        'n_instruments': K,
        'weighted_mean_excess': T_obs,
        'p_two_sided': p_two,
        'p_site2_loss': p_loss,
        'direction': 'site2_loss' if T_obs < 0 else 'site2_gain',
        'frac_instruments_consistent': n_consistent / K,
        'null_mean': m,
        'null_std': float(null_T.std()),
        'per_instrument': per_inst,
    }


def _secondary_free_mask(raw_calls_seq, lo0, hi0, exclude_ids=()):
    """§2.5 control. Boolean mask over reads: True = the read carries
    NO raw (pre-consensus) SNV located within the site-2 window
    [lo0, hi0] (0-based ref), excluding the instrument's own site-1
    call ids. A read with a competing site-2-local mutation could
    explain a distal site-2 change by itself (double-mutant /
    barcode-collision haplotype) rather than cooperativity, so the
    call must survive when those reads are removed.
    """
    excl = set(exclude_ids)
    n = len(raw_calls_seq)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        for vid in raw_calls_seq[i]:
            if vid in excl:
                continue
            pos1, _r, _a, _ct = parse_variant_id_fields(vid)
            if pos1 is None:
                continue
            if lo0 <= (pos1 - 1) <= hi0:    # 1-based id -> 0-based ref
                keep[i] = False
                break
    return keep


def run_cooccupancy(rd, site_pairs, n_iter, random_seed):
    """Driver: for each (site1, site2) pair, build the WT-conditional
    bootstrap ONCE (reused across the site-1's instruments), aggregate
    across instruments, BH across all tested pairs, then mark
    co-occupancy CALLS = BH q < call_q AND directional consistency
    across independent instruments >= min_consistency (the Phase-2
    sign-consistency analogue: cooperativity requires independent SNVs
    to AGREE on direction, not merely a nonzero mean excess). Returns
    pair results sorted by q."""
    cfg = site_pairs['cfg']
    occ_frac = cfg['site_occupancy_frac']
    min_stratum = cfg['min_stratum']
    sites_by_id = site_pairs['sites_by_id']
    instruments = site_pairs['site1_instruments']
    results = []
    for s1_id, s2_id in site_pairs['pairs']:
        s1, s2 = sites_by_id[s1_id], sites_by_id[s2_id]
        vids = [v for v in instruments.get(s1_id, [])
                if v in rd.variant_indices]
        if not vids:
            continue
        rng = np.random.default_rng(
            stable_variant_seed(random_seed, f"cooc_{s1_id}|{s2_id}"))
        wt_cond = wt_conditional_bootstrap(
            rd, s1, s2, rd.wt_indices, occ_frac, n_iter, rng,
            min_stratum)
        if wt_cond is None:
            continue
        inst_states = []
        inst_vi = []
        for vid in vids:
            vi = rd.variant_indices[vid]
            inst_vi.append((vid, vi))
            inst_states.append((
                read_site_occupied(rd, s1, vi, occ_frac),
                read_site_occupied(rd, s2, vi, occ_frac)))
        agg = _cooccupancy_aggregate_from_states(
            inst_states, wt_cond, rng)
        if agg is None:
            continue
        agg.update({
            'site1_id': s1_id, 'site2_id': s2_id,
            'site1_bin': s1['bin'], 'site2_bin': s2['bin'],
            'site1_abs': (s1['abs_start'], s1['abs_end']),
            'site2_abs': (s2['abs_start'], s2['abs_end']),
            'wt_p2_given_1': wt_cond['p1'],
            'wt_p2_given_0': wt_cond['p0'],
            'instrument_variant_ids': vids,
        })

        # ---- §2.5 secondary-mutation control (P3.5) ----
        # Re-test on reads carrying NO competing site-2-local raw SNV;
        # a genuine cooperative call survives this, a double-mutant /
        # haplotype artifact collapses. Cross-instrument: report the
        # max frequency of any single such secondary SNV.
        if cfg['secondary_screen']:
            W = cfg['secondary_window']
            lo0 = s2['abs_start'] - W
            hi0 = s2['abs_end'] + W
            sec_counts = Counter()
            n_pool = 0
            sf_states = []
            for (vid, vi), (o1, o2) in zip(inst_vi, inst_states):
                raw_seq = rd.raw_calls[vi]
                keep = _secondary_free_mask(raw_seq, lo0, hi0,
                                            exclude_ids=(vid,))
                n_pool += len(vi)
                for i in range(len(vi)):
                    for rv in raw_seq[i]:
                        if rv == vid:
                            continue
                        p1pos, _r, _a, _c = parse_variant_id_fields(rv)
                        if p1pos is not None and lo0 <= p1pos - 1 <= hi0:
                            sec_counts[rv] += 1
                if keep.sum() >= cfg['min_stratum']:
                    sf_states.append((np.asarray(o1)[keep],
                                      np.asarray(o2)[keep]))
            max_sec_freq = (max(sec_counts.values()) / n_pool
                            if sec_counts and n_pool else 0.0)
            sf_rng = np.random.default_rng(stable_variant_seed(
                random_seed, f"cooc_secfree_{s1_id}|{s2_id}"))
            sf_agg = (_cooccupancy_aggregate_from_states(
                          sf_states, wt_cond, sf_rng)
                      if len(sf_states) >= cfg['min_variant_instruments']
                      else None)
            agg['secondary_screen'] = {
                'max_secondary_freq': float(max_sec_freq),
                'n_secondary_distinct': len(sec_counts),
                'n_instruments_secfree': len(sf_states),
                'reads_total': int(n_pool),
                'secfree_excess': (sf_agg['weighted_mean_excess']
                                   if sf_agg else None),
                'secfree_p_two_sided': (sf_agg['p_two_sided']
                                        if sf_agg else None),
                'secfree_frac_consistent': (
                    sf_agg['frac_instruments_consistent']
                    if sf_agg else None),
            }
        results.append(agg)
    if results:
        qs = benjamini_hochberg(np.array(
            [r['p_two_sided'] for r in results]))
        cq, mc = cfg['call_q'], cfg['min_consistency']
        for r, q in zip(results, qs):
            r['fdr_q'] = float(q)
            r['is_call'] = bool(
                q < cq
                and r['frac_instruments_consistent'] >= mc)
            # §2.5: a call is VALIDATED only if it survives the
            # secondary-mutation screen — the secondary-free re-test is
            # still significant, consistent, and SAME direction.
            ss = r.get('secondary_screen')
            if not cfg['secondary_screen']:
                survives = True
            elif not r['is_call']:
                survives = False
            elif ss is None or ss['secfree_p_two_sided'] is None:
                survives = False        # could not confirm -> not valid
            else:
                survives = bool(
                    ss['secfree_p_two_sided'] < cq
                    and ss['secfree_frac_consistent'] >= mc
                    and (np.sign(ss['secfree_excess'])
                         == np.sign(r['weighted_mean_excess'])))
            r['secondary_artifact'] = bool(r['is_call']
                                           and not survives)
            r['validated_call'] = bool(r['is_call'] and survives)
        results.sort(key=lambda r: r['p_two_sided'])
    n_call = sum(1 for r in results if r.get('is_call'))
    n_valid = sum(1 for r in results if r.get('validated_call'))
    n_artifact = sum(1 for r in results if r.get('secondary_artifact'))
    logging.info(
        f"  Co-occupancy P3.4/3.5: tested {len(results)} "
        f"(site1,site2) pairs; FDR<{cfg['call_q']}="
        f"{sum(1 for r in results if r.get('fdr_q', 1) < cfg['call_q'])}"
        f"; calls={n_call}; VALIDATED (survive §2.5 secondary "
        f"screen)={n_valid}; secondary-mutation artifacts="
        f"{n_artifact}")
    return results
