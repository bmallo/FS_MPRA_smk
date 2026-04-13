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
import json
import logging
import re
import struct
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import shared_memory

import numpy as np
import pysam
from scipy import ndimage

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

CLUSTER_WIDTH_RANGES = {
    'sub-TF_10-19bp':      (5, 30),
    'TF_20-40bp':          (10, 60),
    'PIC_41-80bp':         (25, 100),
    'Nucleosome_81plusbp':  (50, 200),
}

TAG_NUC_COUNT = 'nc'
TAG_PROMOTER_VARIANT = 'PV'
TAG_VARIANT_COUNT = 'VC'
TAG_FP_STARTS = 'ns'
TAG_FP_LENGTHS = 'nl'
TAG_FP_QUAL = 'nq'


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
            self._coverage_rows = {label: [] for label in bin_labels}

        # Populated on freeze()
        self.nuc_counts = None
        self.variant_ids = None
        self.variant_counts = None
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

        for label in self.bin_labels:
            new_mat = np.zeros((new_cap, self.analysis_length), dtype=np.uint8)
            new_mat[:self._capacity] = self._coverage_matrices[label]
            self._coverage_matrices[label] = new_mat
        self._capacity = new_cap

    def add_read_coverage(self, nuc_count, variant_id, variant_count,
                          bin_coverage):
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
        for label, cov in bin_coverage.items():
            if label in self.bin_to_idx:
                self._coverage_matrices[label][i] = cov
        # Bins not in bin_coverage stay as zeros (initialized in __init__)
        self.n_reads += 1

    def add_read(self, nuc_count, variant_id, variant_count, footprints,
                 ref_length):
        """Legacy add_read: append-based, used when analysis_length not set."""
        assert not self._preallocated and not self._frozen
        self._nuc_counts.append(nuc_count if nuc_count is not None else -1)
        self._variant_ids.append(variant_id)
        self._variant_counts.append(variant_count)
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


def parse_bam(bam_path, size_bins, bin_labels, target_chrom=None,
              analysis_region=None, estimated_reads=None):
    """Parse a tagged BAM into a ReadData object.

    Parameters
    ----------
    bam_path : str
        Path to BAM file (must have PV/VC tags from stage 2).
    size_bins : list of (min, max, label)
        Footprint size bin definitions.
    bin_labels : list of str
        Bin label names.
    target_chrom : str or None
        Target chromosome name (auto-detected if single-ref BAM).
    analysis_region : tuple of (start, end) or None
        If provided, coverage matrices are restricted to this reference
        region, using the optimized CIGAR-walk q2r and vectorized
        footprint conversion. Reduces memory ~15x for typical cases.
        If None, falls back to legacy full-width behavior.
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
    bin_tuples = [(b[0], b[1], b[2]) for b in size_bins]

    # Pre-compute bin boundaries for optimized path
    if analysis_region is not None:
        bin_boundaries = [(b[0], b[1]) for b in size_bins]

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        ref_info = {sq['SN']: sq['LN'] for sq in bam.header['SQ']}
        logging.info(f"References in BAM: {ref_info}")
        ref_name, ref_length = resolve_target_chrom(ref_info, target_chrom)
        logging.info(f"Using reference: {ref_name} ({ref_length:,} bp)")

        # Estimate read count for pre-allocation
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

        # Choose construction mode
        if analysis_region is not None:
            a_start, a_end = analysis_region
            a_length = a_end - a_start
            logging.info(f"Analysis region: {a_start}-{a_end} "
                         f"({a_length} bp of {ref_length} bp)")
            rd = ReadData(ref_length, bin_labels,
                          analysis_length=a_length,
                          analysis_start=a_start,
                          estimated_reads=estimated_reads)
        else:
            rd = ReadData(ref_length, bin_labels)

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
            try:
                fp_starts = read.get_tag(TAG_FP_STARTS)
                fp_lengths = read.get_tag(TAG_FP_LENGTHS)
            except KeyError:
                stats['missing_fp_tags'] += 1
                continue

            if pv_value == "WT" or vc_value == 0:
                variant_id = "WT"
            else:
                variants = parse_variant_tag(pv_value)
                variant_id = (variants[0] if len(variants) == 1
                              else json.dumps(variants))

            if analysis_region is not None:
                # Optimized path: CIGAR-walk q2r + vectorized coverage
                if len(fp_starts) == 0:
                    stats['no_footprints'] += 1
                q2r_arr = build_q2r_array(read)
                bin_cov = convert_footprints_to_coverage(
                    fp_starts, fp_lengths, q2r_arr,
                    a_start, a_length, bin_boundaries, bin_labels)
                rd.add_read_coverage(nuc_count, variant_id, int(vc_value),
                                     bin_cov)
            else:
                # Legacy path
                fp_starts = parse_tag_array(fp_starts)
                fp_lengths = parse_tag_array(fp_lengths)
                if len(fp_starts) == 0:
                    stats['no_footprints'] += 1
                fp_quals = None
                try:
                    fp_quals = parse_tag_array(read.get_tag(TAG_FP_QUAL))
                except KeyError:
                    pass
                q2r = build_query_to_ref_map(
                    read.get_aligned_pairs(matches_only=False))
                footprints = convert_footprints(
                    fp_starts, fp_lengths, fp_quals, q2r,
                    ref_length, bin_tuples)
                rd.add_read(nuc_count, variant_id, int(vc_value),
                            footprints, ref_length)
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

def _detect_clusters_adaptive(abs_delta, threshold, gap_tolerance=2,
                              merge_distance=5, min_width=3):
    sig = abs_delta > threshold
    n = len(sig)
    if not np.any(sig):
        return []
    if gap_tolerance > 0:
        bridged = sig.copy()
        for i in range(n):
            if not sig[i]:
                left = np.any(sig[max(0, i - gap_tolerance):i])
                right = np.any(sig[i + 1:min(n, i + gap_tolerance + 1)])
                if left and right:
                    bridged[i] = True
        sig = bridged
    labeled, n_clusters = ndimage.label(sig)
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
        rd = abs_delta[s:e + 1]
        clusters.append({
            'start': s, 'end': e, 'width': w,
            'sum_abs_delta': float(np.sum(rd)),
            'max_abs_delta': float(np.max(rd)),
            'mean_abs_delta': float(np.mean(rd)),
            'peak_position': s + int(np.argmax(rd)),
        })
    return clusters


def _detect_clusters_from_mask(sig_mask, abs_delta, signed_delta,
                               gap_tolerance=2, merge_distance=5,
                               min_width=3):
    n = len(sig_mask)
    if not np.any(sig_mask):
        return []
    if gap_tolerance > 0:
        bridged = sig_mask.copy()
        for i in range(n):
            if not sig_mask[i]:
                left = np.any(sig_mask[max(0, i - gap_tolerance):i])
                right = np.any(sig_mask[i + 1:min(n, i + gap_tolerance + 1)])
                if left and right:
                    bridged[i] = True
        sig_mask = bridged
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
        rs = signed_delta[s:e + 1]
        clusters.append({
            'start': s, 'end': e, 'width': w,
            'sum_abs_delta': float(np.sum(ra)),
            'max_abs_delta': float(np.max(ra)),
            'mean_abs_delta': float(np.mean(ra)),
            'mean_signed_delta': float(np.mean(rs)),
            'peak_position': s + int(np.argmax(ra)),
            'direction': 'loss' if np.mean(rs) < 0 else 'gain',
        })
    return clusters


# ============================================================================
# Null Calibration — Worker
# ============================================================================

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
    cluster_stats = []

    for li, seed in enumerate(it_seeds):
        rng = np.random.default_rng(seed)
        pv_idx = _nc_matched_batch_draw(
            nc_draw_counts, wt_nc_pools, 1, subsample_size, rng)[0]
        iter_cl = {}
        for bi, label in enumerate(bin_labels):
            pv_occ = shared_mats[label][pv_idx].mean(axis=0)
            delta = pv_occ - wt_occ[label]
            delta_out[li, bi, :] = delta.astype(np.float32)
            abs_d = np.abs(delta)
            thresh = 0.0
            if cluster_threshold_quantile is not None:
                thresh = np.percentile(abs_d, cluster_threshold_quantile * 100)
            if absolute_delta_threshold is not None:
                thresh = max(thresh, absolute_delta_threshold)
            cl = _detect_clusters_adaptive(abs_d, thresh, gap_tolerance,
                                           merge_distance)
            iter_cl[label] = {
                'cluster_sums': [c['sum_abs_delta'] for c in cl],
                'max_cluster_sum': max((c['sum_abs_delta'] for c in cl),
                                       default=0.0),
            }
        cluster_stats.append(iter_cl)

    for shm in shm_handles:
        shm.close()
    return delta_out, cluster_stats


# ============================================================================
# Null Calibration — Main
# ============================================================================

def run_null_calibration(rd, wt_idx, ground_truth_nc, coverage_level,
                         min_nuc=None, max_nuc=None,
                         n_iterations=2000, random_seed=42,
                         analysis_region=None,
                         cluster_threshold_quantile=0.95,
                         absolute_delta_threshold=None,
                         gap_tolerance=2, merge_distance=5,
                         n_workers=1):
    """Generate empirical null at a given coverage level."""
    logging.info(f"  Null calibration: N={coverage_level}, "
                 f"iters={n_iterations}, workers={n_workers}")

    bin_labels = rd.bin_labels
    n_bins = len(bin_labels)
    if analysis_region is not None:
        a_start, a_end = analysis_region
    else:
        a_start, a_end = 0, rd.ref_length
    analysis_length = a_end - a_start

    wt_f = rd.get_indices_filtered(wt_idx, min_nuc=min_nuc, max_nuc=max_nuc)
    # If matrices are already analysis-width, row-slice only; otherwise
    # also slice columns to the analysis region.
    needs_col_slice = (rd.analysis_length != analysis_length
                       or rd.analysis_start != a_start)
    if needs_col_slice:
        col_indices = np.arange(a_start, a_end)
        wt_mats = {l: rd.coverage_matrices[l][np.ix_(wt_f, col_indices)]
                   for l in bin_labels}
    else:
        wt_mats = {l: rd.coverage_matrices[l][wt_f]
                   for l in bin_labels}
    wt_occ = {l: wt_mats[l].mean(axis=0).astype(np.float64)
              for l in bin_labels}

    wt_nc = rd.nuc_counts[wt_f]
    nc_vals = ground_truth_nc['nc_vals']
    nc_fracs = ground_truth_nc['nc_fracs']
    wt_nc_pools = {}
    for nc_val in nc_vals:
        pool = np.where(wt_nc == nc_val)[0]
        if len(pool) > 0:
            wt_nc_pools[nc_val] = pool
    wt_nc_pools_ser = list(wt_nc_pools.items())

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
    all_cluster_stats = []

    if n_workers > 1 and n_iterations > 1:
        shm_info, shm_objects = create_shared_matrices(wt_mats, bin_labels)
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
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = {executor.submit(_null_worker_shared, wa): ci
                           for ci, wa in enumerate(worker_args)}
                results_by_chunk = {}
                for future in as_completed(futures):
                    results_by_chunk[futures[future]] = future.result()
                offset = 0
                for ci in range(len(seed_chunks)):
                    dc, cc = results_by_chunk[ci]
                    n_c = dc.shape[0]
                    null_delta[offset:offset + n_c] = dc
                    all_cluster_stats.extend(cc)
                    offset += n_c
        finally:
            cleanup_shared_memory(shm_objects)
    else:
        for it in tqdm(range(n_iterations), desc="    Null iters",
                       disable=not HAS_TQDM):
            rng = np.random.default_rng(all_seeds[it])
            pv_idx = _nc_matched_batch_draw(
                nc_draw_counts, wt_nc_pools, 1, coverage_level, rng)[0]
            iter_cl = {}
            for bi, label in enumerate(bin_labels):
                pv_occ = wt_mats[label][pv_idx].mean(axis=0)
                delta = pv_occ - wt_occ[label]
                null_delta[it, bi, :] = delta.astype(np.float32)
                abs_d = np.abs(delta)
                thresh = 0.0
                if cluster_threshold_quantile is not None:
                    thresh = np.percentile(abs_d,
                                           cluster_threshold_quantile * 100)
                if absolute_delta_threshold is not None:
                    thresh = max(thresh, absolute_delta_threshold)
                cl = _detect_clusters_adaptive(abs_d, thresh, gap_tolerance,
                                               merge_distance)
                iter_cl[label] = {
                    'cluster_sums': [c['sum_abs_delta'] for c in cl],
                    'max_cluster_sum': max((c['sum_abs_delta'] for c in cl),
                                           default=0.0),
                }
            all_cluster_stats.append(iter_cl)

    # Summarize
    results = {
        'null_delta': null_delta, 'wt_occ': wt_occ,
        'coverage_level': coverage_level, 'n_iterations': n_iterations,
        'analysis_region': (a_start, a_end),
        'analysis_length': analysis_length,
        'summary': {}, 'null_cluster_sums': {},
        'null_max_cluster_sums': {},
    }
    for bi, label in enumerate(bin_labels):
        d = null_delta[:, bi, :]
        results['summary'][label] = {
            'pos_null_mean': d.mean(axis=0),
            'pos_null_std': d.std(axis=0),
        }
        sums = []
        max_sums = []
        for cs in all_cluster_stats:
            if label in cs:
                sums.extend(cs[label]['cluster_sums'])
                max_sums.append(cs[label]['max_cluster_sum'])
        results['null_cluster_sums'][label] = np.array(sums, dtype=np.float32)
        results['null_max_cluster_sums'][label] = np.array(max_sums,
                                                            dtype=np.float32)
        p95 = float(np.percentile(np.abs(d).flatten(), 95))
        logging.info(f"    {label}: null |delta| p95={p95:.4f}")

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

def test_single_variant(rd, wt_idx, var_indices, var_id, ground_truth_nc,
                        null_results_by_depth, analysis_region,
                        promoter_start, promoter_end,
                        min_nuc=None, max_nuc=None, random_seed=42,
                        cluster_threshold_quantile=0.95,
                        absolute_delta_threshold=None,
                        gap_tolerance=2, merge_distance=5):
    """Test one variant. Returns a dict of per-bin results."""
    bin_labels = rd.bin_labels
    a_start, a_end = analysis_region
    analysis_length = a_end - a_start
    prom_rel_s = max(0, promoter_start - a_start)
    prom_rel_e = min(analysis_length, promoter_end - a_start)
    prom_slice = slice(prom_rel_s, prom_rel_e)

    nc_vals = ground_truth_nc['nc_vals']
    nc_fracs = ground_truth_nc['nc_fracs']
    rng = np.random.default_rng(random_seed)

    var_matched, var_n = nc_matched_subsample(
        var_indices, rd.nuc_counts, nc_fracs, nc_vals, rng, replace=False)

    if var_n < 20:
        return None

    available_depths = sorted(null_results_by_depth.keys())
    closest = min(available_depths, key=lambda d: abs(d - var_n))
    null_res = null_results_by_depth[closest]

    wt_occ = null_res['wt_occ']
    # If matrices are already analysis-width, row-slice only
    needs_col_slice = (rd.analysis_length != analysis_length
                       or rd.analysis_start != a_start)
    if needs_col_slice:
        col_indices = np.arange(a_start, a_end)

    results = {'variant_id': var_id, 'n_raw': len(var_indices),
               'n_nc_matched': var_n, 'null_depth_used': closest}

    best_cluster_p = 1.0

    for bi, label in enumerate(bin_labels):
        if needs_col_slice:
            var_occ = rd.coverage_matrices[label][
                np.ix_(var_matched, col_indices)].mean(axis=0)
        else:
            var_occ = rd.coverage_matrices[label][
                var_matched].mean(axis=0)
        delta_obs = var_occ - wt_occ[label]

        null_delta = null_res['null_delta'][:, bi, :]
        n_null = null_delta.shape[0]
        abs_obs = np.abs(delta_obs)
        abs_null = np.abs(null_delta)

        exceed = np.sum(abs_null >= abs_obs[np.newaxis, :], axis=0)
        empirical_p = (exceed + 1) / (n_null + 1)

        prom_p = empirical_p[prom_slice]
        prom_q = benjamini_hochberg(prom_p)
        q_values = np.ones(analysis_length)
        q_values[prom_slice] = prom_q

        null_mean = null_res['summary'][label]['pos_null_mean']
        null_std = np.maximum(null_res['summary'][label]['pos_null_std'], 1e-6)
        z_scores = (delta_obs - null_mean) / null_std

        pos_thresh = np.percentile(abs_null,
                                   cluster_threshold_quantile * 100, axis=0)
        sig_mask = abs_obs > pos_thresh
        if absolute_delta_threshold is not None:
            sig_mask &= abs_obs > absolute_delta_threshold

        clusters = _detect_clusters_from_mask(
            sig_mask, abs_obs, delta_obs, gap_tolerance, merge_distance)

        null_csums = null_res['null_cluster_sums'].get(label, np.array([]))
        null_max_csums = null_res['null_max_cluster_sums'].get(
            label, np.array([]))

        sig_clusters = []
        for c in clusters:
            if c['end'] < prom_rel_s or c['start'] >= prom_rel_e:
                continue
            if len(null_csums) > 0:
                c['sum_p'] = max(
                    float(np.mean(null_csums >= c['sum_abs_delta'])),
                    1.0 / (len(null_csums) + 1))
            else:
                c['sum_p'] = np.nan
            if len(null_max_csums) > 0:
                c['max_sum_p'] = float(
                    np.mean(null_max_csums >= c['sum_abs_delta']))
            else:
                c['max_sum_p'] = np.nan
            c['abs_start'] = c['start'] + a_start
            c['abs_end'] = c['end'] + a_start
            if c.get('sum_p', 1.0) < 0.05:
                sig_clusters.append(c)
            if c.get('sum_p', 1.0) < best_cluster_p:
                best_cluster_p = c.get('sum_p', 1.0)

        n_sig_pos = int(np.sum(q_values[prom_slice] < 0.10))

        results[label] = {
            'delta_obs': delta_obs.astype(np.float32),
            'variant_occ': var_occ.astype(np.float32),
            'empirical_p': empirical_p.astype(np.float32),
            'q_values': q_values.astype(np.float32),
            'z_scores': z_scores.astype(np.float32),
            'n_sig_positions_fdr10': n_sig_pos,
            'max_abs_delta': float(np.max(abs_obs)),
            'significant_clusters': sig_clusters,
            'all_promoter_clusters': [c for c in clusters
                                      if c['end'] >= prom_rel_s
                                      and c['start'] < prom_rel_e],
        }

    results['best_cluster_p'] = best_cluster_p
    return results


# ============================================================================
# Parallel Per-Variant Testing
# ============================================================================

def _variant_test_worker(args):
    """Worker for parallel variant testing via shared memory."""
    (variant_batch, shm_cov_info, shm_nc_info, bin_labels,
     null_results_by_depth_ser, analysis_region,
     promoter_start, promoter_end, ground_truth_nc,
     min_nuc, max_nuc, random_seed_base,
     cluster_threshold_quantile, absolute_delta_threshold,
     gap_tolerance, merge_distance, analysis_length, analysis_start) = args

    # Attach to shared memory for coverage matrices and nuc_counts
    cov_mats = {}
    shm_handles = []
    for label in bin_labels:
        arr, shm = load_shared_matrix(shm_cov_info[label])
        cov_mats[label] = arr
        shm_handles.append(shm)

    nc_arr, nc_shm = load_shared_matrix(shm_nc_info)
    shm_handles.append(nc_shm)

    a_start, a_end = analysis_region
    prom_rel_s = max(0, promoter_start - a_start)
    prom_rel_e = min(analysis_length, promoter_end - a_start)
    prom_slice = slice(prom_rel_s, prom_rel_e)

    nc_vals = ground_truth_nc['nc_vals']
    nc_fracs = ground_truth_nc['nc_fracs']

    # Deserialize null results (just the summaries, not full delta arrays)
    null_results_by_depth = null_results_by_depth_ser

    batch_results = []
    for vid, var_indices, seed in variant_batch:
        rng = np.random.default_rng(seed)
        var_matched, var_n = nc_matched_subsample(
            var_indices, nc_arr, nc_fracs, nc_vals, rng, replace=False)

        if var_n < 20:
            batch_results.append(None)
            continue

        available_depths = sorted(null_results_by_depth.keys())
        closest = min(available_depths, key=lambda d: abs(d - var_n))
        null_res = null_results_by_depth[closest]
        wt_occ = null_res['wt_occ']

        results = {'variant_id': vid, 'n_raw': len(var_indices),
                   'n_nc_matched': var_n, 'null_depth_used': closest}
        best_cluster_p = 1.0

        for bi, label in enumerate(bin_labels):
            var_occ = cov_mats[label][var_matched].mean(axis=0)
            delta_obs = var_occ - wt_occ[label]

            null_delta = null_res['null_delta'][:, bi, :]
            n_null = null_delta.shape[0]
            abs_obs = np.abs(delta_obs)
            abs_null = np.abs(null_delta)

            exceed = np.sum(abs_null >= abs_obs[np.newaxis, :], axis=0)
            empirical_p = (exceed + 1) / (n_null + 1)

            prom_p = empirical_p[prom_slice]
            prom_q = benjamini_hochberg(prom_p)
            q_values = np.ones(analysis_length)
            q_values[prom_slice] = prom_q

            null_mean = null_res['summary'][label]['pos_null_mean']
            null_std = np.maximum(
                null_res['summary'][label]['pos_null_std'], 1e-6)
            z_scores = (delta_obs - null_mean) / null_std

            pos_thresh = np.percentile(
                abs_null, cluster_threshold_quantile * 100, axis=0)
            sig_mask = abs_obs > pos_thresh
            if absolute_delta_threshold is not None:
                sig_mask &= abs_obs > absolute_delta_threshold

            clusters = _detect_clusters_from_mask(
                sig_mask, abs_obs, delta_obs, gap_tolerance, merge_distance)

            null_csums = null_res['null_cluster_sums'].get(
                label, np.array([]))
            null_max_csums = null_res['null_max_cluster_sums'].get(
                label, np.array([]))

            sig_clusters = []
            for c in clusters:
                if c['end'] < prom_rel_s or c['start'] >= prom_rel_e:
                    continue
                if len(null_csums) > 0:
                    c['sum_p'] = max(
                        float(np.mean(null_csums >= c['sum_abs_delta'])),
                        1.0 / (len(null_csums) + 1))
                else:
                    c['sum_p'] = np.nan
                if len(null_max_csums) > 0:
                    c['max_sum_p'] = float(
                        np.mean(null_max_csums >= c['sum_abs_delta']))
                else:
                    c['max_sum_p'] = np.nan
                c['abs_start'] = c['start'] + a_start
                c['abs_end'] = c['end'] + a_start
                if c.get('sum_p', 1.0) < 0.05:
                    sig_clusters.append(c)
                if c.get('sum_p', 1.0) < best_cluster_p:
                    best_cluster_p = c.get('sum_p', 1.0)

            n_sig_pos = int(np.sum(q_values[prom_slice] < 0.10))

            results[label] = {
                'delta_obs': delta_obs.astype(np.float32),
                'variant_occ': var_occ.astype(np.float32),
                'empirical_p': empirical_p.astype(np.float32),
                'q_values': q_values.astype(np.float32),
                'z_scores': z_scores.astype(np.float32),
                'n_sig_positions_fdr10': n_sig_pos,
                'max_abs_delta': float(np.max(abs_obs)),
                'significant_clusters': sig_clusters,
                'all_promoter_clusters': [c for c in clusters
                                          if c['end'] >= prom_rel_s
                                          and c['start'] < prom_rel_e],
            }

        results['best_cluster_p'] = best_cluster_p
        batch_results.append(results)

    for shm in shm_handles:
        shm.close()
    return batch_results


def run_variant_testing_parallel(rd, wt_idx, variant_groups, ground_truth_nc,
                                 null_results_by_depth, analysis_region,
                                 promoter_start, promoter_end,
                                 min_nuc=None, max_nuc=None,
                                 random_seed=42,
                                 cluster_threshold_quantile=0.95,
                                 absolute_delta_threshold=None,
                                 gap_tolerance=2, merge_distance=5,
                                 n_workers=1):
    """Test all variants, optionally in parallel using shared memory.

    Returns list of result dicts (same format as test_single_variant).
    """
    a_start, a_end = analysis_region
    analysis_length = a_end - a_start
    n_total = len(variant_groups)

    if n_workers <= 1 or n_total <= 1:
        # Serial path
        all_results = []
        for vid, var_indices in tqdm(variant_groups.items(),
                                     desc="Testing variants",
                                     disable=not HAS_TQDM):
            vr = test_single_variant(
                rd, wt_idx, var_indices, vid, ground_truth_nc,
                null_results_by_depth, analysis_region,
                promoter_start, promoter_end,
                min_nuc=min_nuc, max_nuc=max_nuc,
                random_seed=random_seed + hash(vid) % (2**31),
                cluster_threshold_quantile=cluster_threshold_quantile,
                absolute_delta_threshold=absolute_delta_threshold,
                gap_tolerance=gap_tolerance,
                merge_distance=merge_distance)
            if vr is not None:
                all_results.append(vr)
        return all_results

    # Parallel path: put coverage matrices and nuc_counts in shared memory
    logging.info(f"  Parallel variant testing: {n_total} variants, "
                 f"{n_workers} workers")

    shm_cov_info, shm_cov_objects = create_shared_matrices(
        rd.coverage_matrices, rd.bin_labels)

    # nuc_counts into shared memory
    nc_shm = shared_memory.SharedMemory(
        create=True, size=rd.nuc_counts.nbytes)
    nc_shared = np.ndarray(rd.nuc_counts.shape, dtype=rd.nuc_counts.dtype,
                           buffer=nc_shm.buf)
    nc_shared[:] = rd.nuc_counts[:]
    shm_nc_info = {'name': nc_shm.name, 'shape': rd.nuc_counts.shape,
                   'dtype': str(rd.nuc_counts.dtype)}

    # Build task list: (vid, var_indices, seed)
    tasks = [(vid, var_indices,
              random_seed + hash(vid) % (2**31))
             for vid, var_indices in variant_groups.items()]

    # Chunk into batches
    batch_size = max(1, len(tasks) // (n_workers * 4))
    batches = [tasks[i:i + batch_size]
               for i in range(0, len(tasks), batch_size)]

    worker_args = [
        (batch, shm_cov_info, shm_nc_info, rd.bin_labels,
         null_results_by_depth, analysis_region,
         promoter_start, promoter_end, ground_truth_nc,
         min_nuc, max_nuc, random_seed,
         cluster_threshold_quantile, absolute_delta_threshold,
         gap_tolerance, merge_distance, analysis_length, rd.analysis_start)
        for batch in batches
    ]

    all_results = []
    try:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_variant_test_worker, wa): ci
                       for ci, wa in enumerate(worker_args)}
            results_by_chunk = {}
            for future in as_completed(futures):
                results_by_chunk[futures[future]] = future.result()
            for ci in range(len(batches)):
                for vr in results_by_chunk[ci]:
                    if vr is not None:
                        all_results.append(vr)
    finally:
        cleanup_shared_memory(shm_cov_objects)
        nc_shm.close()
        nc_shm.unlink()

    logging.info(f"  Tested {len(all_results)} / {n_total} variants")
    return all_results
