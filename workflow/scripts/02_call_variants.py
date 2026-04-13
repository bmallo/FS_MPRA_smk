#!/usr/bin/env python3
"""
MPRA Variant-Barcode Tagger

Tags PacBio MPRA reads with promoter variants and barcode sequences.
Includes barcode clustering by edit distance (neighborhood hashing) and
position-level consensus variant calling, inspired by PacRAT/Pacybara.
Generates statistics, diagnostic plots, and QC reports.
"""

import argparse
import pysam
import json
import logging
import sys
import base64
import io
from pathlib import Path
from dataclasses import dataclass, field
from collections import defaultdict, Counter
from typing import Optional
from enum import Enum

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

BASES = 'ACGT'


class ExclusionReason(Enum):
    PASS = "pass"
    UNMAPPED = "unmapped"
    SECONDARY = "secondary_alignment"
    SUPPLEMENTARY = "supplementary_alignment"
    NO_PROMOTER_COVERAGE = "no_promoter_coverage"
    NO_BARCODE_COVERAGE = "no_barcode_coverage"
    BARCODE_TOO_SHORT = "barcode_too_short"
    BARCODE_TOO_LONG = "barcode_too_long"
    LOW_VARIANT_QUALITY = "low_variant_quality"
    FORBIDDEN_CIGAR = "forbidden_cigar_operation"
    HAS_INDEL_VARIANT = "has_indel_variant"


class VariantType(Enum):
    SNV = "snv"
    INSERTION = "insertion"
    DELETION = "deletion"


@dataclass
class VariantCall:
    position: int
    ref: str
    alt: str
    var_type: VariantType
    quality: Optional[int] = None

    @property
    def id(self) -> str:
        if self.var_type == VariantType.SNV:
            return f"{self.position}:{self.ref}>{self.alt}"
        elif self.var_type == VariantType.INSERTION:
            return f"{self.position}:+{self.alt}"
        else:
            return f"{self.position}:{len(self.ref)}{self.ref}"

    @property
    def is_transition(self) -> bool:
        if self.var_type != VariantType.SNV:
            return False
        transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
        return (self.ref.upper(), self.alt.upper()) in transitions


@dataclass
class ProcessingStats:
    total_reads: int = 0
    passed_reads: int = 0
    exclusion_counts: Counter = field(default_factory=Counter)
    wt_reads: int = 0
    single_variant_reads: int = 0
    multi_variant_reads: int = 0
    variant_counts: Counter = field(default_factory=Counter)
    variant_type_counts: Counter = field(default_factory=Counter)
    transition_count: int = 0
    transversion_count: int = 0
    barcode_variant_map: dict = field(default_factory=lambda: defaultdict(Counter))
    barcode_lengths: Counter = field(default_factory=Counter)
    variants_per_read_hist: Counter = field(default_factory=Counter)
    variant_nucleosome_counts: dict = field(default_factory=lambda: defaultdict(Counter))
    snv_position_counts: Counter = field(default_factory=Counter)
    ins_position_counts: Counter = field(default_factory=Counter)
    del_position_counts: Counter = field(default_factory=Counter)
    observed_snvs: set = field(default_factory=set)
    observed_insertions: set = field(default_factory=set)
    observed_deletions: set = field(default_factory=set)
    single_variant_snvs: set = field(default_factory=set)


@dataclass
class BarcodeClusterResolution:
    consensus_variant_tag: str
    cluster_barcode: str
    cluster_size: int
    total_reads: int
    is_ambiguous: bool
    variant_support: float
    is_low_confidence: bool
    member_barcodes: list = field(default_factory=list)
    position_details: dict = field(default_factory=dict)


@dataclass
class ClusteringStats:
    total_barcodes_before: int = 0
    total_clusters_after: int = 0
    barcodes_merged: int = 0
    singleton_clusters: int = 0
    cluster_size_distribution: Counter = field(default_factory=Counter)
    resolved_clusters: int = 0
    ambiguous_clusters: int = 0
    low_confidence_clusters: int = 0


@dataclass
class ClusterCorrection:
    """Per-cluster summary of corrections applied during Pass 2."""
    centroid: str
    cluster_size: int
    total_reads: int
    raw_variant_tag: str  # most common raw tag in the cluster
    consensus_variant_tag: str
    reads_corrected: int
    correction_rate: float
    variants_removed: list = field(default_factory=list)   # variant IDs dropped from consensus
    variants_added: list = field(default_factory=list)     # variant IDs in consensus but not in some raw calls
    removed_variant_freqs: dict = field(default_factory=dict)  # var_id -> fraction of reads that had it


@dataclass
class ExpectedCoverage:
    promoter_length: int
    expected_snvs: int
    observed_snvs: int
    snv_coverage: float
    expected_insertions: int
    observed_insertions: int
    insertion_coverage: float
    expected_deletions: int
    observed_deletions: int
    deletion_coverage: float


# ---------------------------------------------------------------------------
# Variant calling from CIGAR
# ---------------------------------------------------------------------------

def parse_promoter_variants(read, ref_seq, region_start_0based, region_end_0based, min_quality=0):
    """Parse CIGAR string to extract variants in the promoter region."""
    variants = []
    if read.reference_end is None or read.reference_start > region_end_0based or read.reference_end <= region_start_0based:
        return [], ExclusionReason.NO_PROMOTER_COVERAGE
    cigar = read.cigartuples
    if cigar is None:
        return [], ExclusionReason.NO_PROMOTER_COVERAGE
    query_seq = read.query_sequence
    query_quals = read.query_qualities
    query_idx = 0
    ref_pos = read.reference_start

    for op, length in cigar:
        if ref_pos > region_end_0based:
            break
        if op == 0:
            return [], ExclusionReason.FORBIDDEN_CIGAR
        elif op == 7:
            query_idx += length
            ref_pos += length
        elif op == 8:
            if ref_pos + length > region_start_0based and ref_pos <= region_end_0based:
                overlap_start = max(ref_pos, region_start_0based)
                overlap_end = min(ref_pos + length, region_end_0based + 1)
                block_offset = overlap_start - ref_pos
                for i in range(overlap_end - overlap_start):
                    abs_ref_pos = overlap_start + i
                    rel_ref_pos = abs_ref_pos - region_start_0based
                    q_idx = query_idx + block_offset + i
                    if rel_ref_pos < len(ref_seq) and q_idx < len(query_seq):
                        ref_base = ref_seq[rel_ref_pos]
                        alt_base = query_seq[q_idx]
                        qual = query_quals[q_idx] if query_quals is not None else None
                        if qual is None or qual >= min_quality:
                            variants.append(VariantCall(position=abs_ref_pos + 1, ref=ref_base, alt=alt_base, var_type=VariantType.SNV, quality=qual))
            query_idx += length
            ref_pos += length
        elif op == 1:
            if region_start_0based <= ref_pos <= region_end_0based:
                inserted_seq = query_seq[query_idx:query_idx + length]
                min_ins_qual = min(query_quals[query_idx:query_idx + length]) if query_quals is not None else None
                if min_ins_qual is None or min_ins_qual >= min_quality:
                    variants.append(VariantCall(position=ref_pos + 1, ref="", alt=inserted_seq, var_type=VariantType.INSERTION, quality=min_ins_qual))
            query_idx += length
        elif op == 2:
            if ref_pos + length > region_start_0based and ref_pos <= region_end_0based:
                overlap_start = max(ref_pos, region_start_0based)
                overlap_end = min(ref_pos + length, region_end_0based + 1)
                if overlap_start < overlap_end:
                    rel_start = overlap_start - region_start_0based
                    rel_end = overlap_end - region_start_0based
                    deleted_seq = ref_seq[rel_start:rel_end]
                    variants.append(VariantCall(position=overlap_start + 1, ref=deleted_seq, alt="", var_type=VariantType.DELETION, quality=None))
            ref_pos += length
        elif op == 4:
            query_idx += length
        elif op == 5:
            pass
        elif op == 3:
            ref_pos += length
    return variants, ExclusionReason.PASS


# ---------------------------------------------------------------------------
# Barcode extraction
# ---------------------------------------------------------------------------

def extract_barcode_seq(read, ref_start_0based, ref_end_0based, expected_length=15, min_length=13, max_length=17):
    """Extract barcode sequence from read at specified reference coordinates.

    Only includes query bases that align to actual reference positions within
    the barcode region. Inserted bases (query bases with no reference position)
    are excluded — for fixed-length synthetic barcodes, insertions are sequencing
    errors and including them corrupts the barcode identity and wastes edit-distance
    budget during clustering.
    """
    if read.is_unmapped:
        return None, None, ExclusionReason.UNMAPPED
    if read.reference_end is None or read.reference_start > ref_end_0based or read.reference_end <= ref_start_0based:
        return None, None, ExclusionReason.NO_BARCODE_COVERAGE
    aligned_pairs = read.get_aligned_pairs(matches_only=False)
    barcode_query_positions = []
    region_started = False
    for query_pos, ref_pos in aligned_pairs:
        if ref_pos is not None and ref_start_0based <= ref_pos < ref_end_0based:
            region_started = True
            if query_pos is not None:
                barcode_query_positions.append(query_pos)
            # query_pos is None here means a deletion in the barcode region;
            # the position is simply skipped (shorter barcode)
        elif region_started and ref_pos is not None and ref_pos >= ref_end_0based:
            break
        elif not region_started and ref_pos is not None and ref_pos >= ref_end_0based:
            return None, None, ExclusionReason.NO_BARCODE_COVERAGE
        # Inserted bases (ref_pos is None) within the region are intentionally
        # excluded — they are sequencing errors for fixed-length barcodes
    if len(barcode_query_positions) < min_length:
        return (None, None, ExclusionReason.NO_BARCODE_COVERAGE) if len(barcode_query_positions) == 0 else (None, None, ExclusionReason.BARCODE_TOO_SHORT)
    if len(barcode_query_positions) > max_length:
        return None, None, ExclusionReason.BARCODE_TOO_LONG
    barcode_seq = ''.join([read.query_sequence[pos] for pos in barcode_query_positions])
    mean_quality = sum(read.query_qualities[pos] for pos in barcode_query_positions) / len(barcode_query_positions) if read.query_qualities is not None else None
    return barcode_seq, mean_quality, ExclusionReason.PASS


# ---------------------------------------------------------------------------
# Combined promoter variant calling + barcode extraction (single CIGAR walk)
# ---------------------------------------------------------------------------

def parse_promoter_and_barcode(read, ref_seq, prom_start_0, prom_end_0,
                               bc_start_0=None, bc_end_0=None,
                               min_quality=0, min_bc_length=13,
                               max_bc_length=17):
    """Parse CIGAR to extract promoter variants and barcode in one pass.

    Replaces separate parse_promoter_variants() + extract_barcode_seq() calls,
    eliminating the expensive get_aligned_pairs() call for barcode extraction.

    Parameters
    ----------
    read : pysam.AlignedSegment
    ref_seq : str
        Reference sequence for the promoter region (prom_start_0-relative).
    prom_start_0, prom_end_0 : int
        0-based promoter region [start, end] (inclusive on both ends,
        matching parse_promoter_variants convention).
    bc_start_0, bc_end_0 : int or None
        0-based barcode region [start, end) (half-open, matching
        extract_barcode_seq convention). None to skip barcode extraction.
    min_quality : int
        Minimum base quality for variant calls.
    min_bc_length, max_bc_length : int
        Barcode length bounds.

    Returns
    -------
    variants : list of VariantCall
    var_status : ExclusionReason
    barcode_seq : str or None
    barcode_qual : float or None
    bc_status : ExclusionReason
    """
    extract_barcode = bc_start_0 is not None

    # Check promoter coverage
    if (read.reference_end is None
            or read.reference_start > prom_end_0
            or read.reference_end <= prom_start_0):
        bc_status = ExclusionReason.NO_BARCODE_COVERAGE if extract_barcode else ExclusionReason.PASS
        return [], ExclusionReason.NO_PROMOTER_COVERAGE, None, None, bc_status

    cigar = read.cigartuples
    if cigar is None:
        bc_status = ExclusionReason.NO_BARCODE_COVERAGE if extract_barcode else ExclusionReason.PASS
        return [], ExclusionReason.NO_PROMOTER_COVERAGE, None, None, bc_status

    query_seq = read.query_sequence
    query_quals = read.query_qualities
    query_idx = 0
    ref_pos = read.reference_start

    variants = []
    barcode_qpos = []  # query positions aligned to barcode region

    # Determine how far we need to walk
    if extract_barcode:
        walk_end = max(prom_end_0, bc_end_0)
    else:
        walk_end = prom_end_0

    for op, length in cigar:
        if ref_pos > walk_end:
            break

        if op == 0:
            # M (alignment match) — forbidden in pbmm2 output
            return [], ExclusionReason.FORBIDDEN_CIGAR, None, None, ExclusionReason.FORBIDDEN_CIGAR

        elif op == 7:
            # = (sequence match): both query and ref advance, no variants
            if extract_barcode:
                # Collect barcode positions within this block
                bc_overlap_s = max(ref_pos, bc_start_0)
                bc_overlap_e = min(ref_pos + length, bc_end_0)
                if bc_overlap_s < bc_overlap_e:
                    offset = bc_overlap_s - ref_pos
                    for j in range(bc_overlap_e - bc_overlap_s):
                        barcode_qpos.append(query_idx + offset + j)
            query_idx += length
            ref_pos += length

        elif op == 8:
            # X (mismatch): both advance, each position is a variant
            # Promoter variant calling
            if ref_pos + length > prom_start_0 and ref_pos <= prom_end_0:
                overlap_start = max(ref_pos, prom_start_0)
                overlap_end = min(ref_pos + length, prom_end_0 + 1)
                block_offset = overlap_start - ref_pos
                for i in range(overlap_end - overlap_start):
                    abs_ref_pos = overlap_start + i
                    rel_ref_pos = abs_ref_pos - prom_start_0
                    q_idx = query_idx + block_offset + i
                    if rel_ref_pos < len(ref_seq) and q_idx < len(query_seq):
                        ref_base = ref_seq[rel_ref_pos]
                        alt_base = query_seq[q_idx]
                        qual = query_quals[q_idx] if query_quals is not None else None
                        if qual is None or qual >= min_quality:
                            variants.append(VariantCall(
                                position=abs_ref_pos + 1, ref=ref_base,
                                alt=alt_base, var_type=VariantType.SNV,
                                quality=qual))
            # Barcode: mismatches still produce aligned positions
            if extract_barcode:
                bc_overlap_s = max(ref_pos, bc_start_0)
                bc_overlap_e = min(ref_pos + length, bc_end_0)
                if bc_overlap_s < bc_overlap_e:
                    offset = bc_overlap_s - ref_pos
                    for j in range(bc_overlap_e - bc_overlap_s):
                        barcode_qpos.append(query_idx + offset + j)
            query_idx += length
            ref_pos += length

        elif op == 1:
            # I (insertion): query advances, ref does not
            if prom_start_0 <= ref_pos <= prom_end_0:
                inserted_seq = query_seq[query_idx:query_idx + length]
                min_ins_qual = (min(query_quals[query_idx:query_idx + length])
                                if query_quals is not None else None)
                if min_ins_qual is None or min_ins_qual >= min_quality:
                    variants.append(VariantCall(
                        position=ref_pos + 1, ref="", alt=inserted_seq,
                        var_type=VariantType.INSERTION, quality=min_ins_qual))
            # Barcode: insertions are excluded (no ref position)
            query_idx += length

        elif op == 2:
            # D (deletion): ref advances, query does not
            if ref_pos + length > prom_start_0 and ref_pos <= prom_end_0:
                overlap_start = max(ref_pos, prom_start_0)
                overlap_end = min(ref_pos + length, prom_end_0 + 1)
                if overlap_start < overlap_end:
                    rel_start = overlap_start - prom_start_0
                    rel_end = overlap_end - prom_start_0
                    deleted_seq = ref_seq[rel_start:rel_end]
                    variants.append(VariantCall(
                        position=overlap_start + 1, ref=deleted_seq, alt="",
                        var_type=VariantType.DELETION, quality=None))
            # Barcode: deletions skip ref positions (shorter barcode)
            ref_pos += length

        elif op == 4:
            # S (soft clip): query advances
            query_idx += length
        elif op == 5:
            # H (hard clip): neither advances
            pass
        elif op == 3:
            # N (ref skip): ref advances
            ref_pos += length

    # Resolve barcode
    barcode_seq = None
    barcode_qual = None
    if extract_barcode:
        if len(barcode_qpos) == 0:
            bc_status = ExclusionReason.NO_BARCODE_COVERAGE
        elif len(barcode_qpos) < min_bc_length:
            bc_status = ExclusionReason.BARCODE_TOO_SHORT
        elif len(barcode_qpos) > max_bc_length:
            bc_status = ExclusionReason.BARCODE_TOO_LONG
        else:
            barcode_seq = ''.join(query_seq[p] for p in barcode_qpos)
            barcode_qual = (sum(query_quals[p] for p in barcode_qpos)
                            / len(barcode_qpos)
                            if query_quals is not None else None)
            bc_status = ExclusionReason.PASS
    else:
        bc_status = ExclusionReason.PASS

    return variants, ExclusionReason.PASS, barcode_seq, barcode_qual, bc_status


# ---------------------------------------------------------------------------
# Barcode clustering by neighborhood hashing
# ---------------------------------------------------------------------------

def generate_levenshtein_neighbors_d1(seq):
    """Generate all sequences within Levenshtein distance 1 of seq."""
    neighbors = set()
    n = len(seq)
    for i in range(n):
        for b in BASES:
            if b != seq[i]:
                neighbors.add(seq[:i] + b + seq[i+1:])
    for i in range(n):
        neighbors.add(seq[:i] + seq[i+1:])
    for i in range(n + 1):
        for b in BASES:
            neighbors.add(seq[:i] + b + seq[i:])
    neighbors.discard(seq)
    return neighbors


def cluster_barcodes(barcode_variant_map, barcode_read_variants,
                     max_edit_distance=2, expected_barcode_length=15,
                     min_jaccard=0.5, quiet=False):
    """
    Cluster barcodes by edit distance using neighborhood hashing,
    with variant-aware merge filtering and length-scaled edit distance.
    
    Three safeguards against spurious merges:
    1. Neighborhood hashing for efficient edit distance lookup
    2. Edit distance threshold scaled by barcode length: barcodes shorter
       than expected get a reduced threshold (max 1 for len < expected)
    3. Variant-aware filtering: a barcode only merges into a cluster if
       the Jaccard index of its variant calls vs the *centroid's* variant
       set exceeds min_jaccard. The centroid's variant set is fixed at
       cluster creation time and never updated, preventing variant set
       drift as members are added.
    
    Jaccard index for variant sets:
      - WT (empty set) only matches WT
      - Two non-empty sets: |intersection| / |union| >= min_jaccard
    """
    barcode_reads = {}
    for bc, var_counts in barcode_variant_map.items():
        barcode_reads[bc] = sum(var_counts.values())

    sorted_barcodes = sorted(barcode_reads.keys(), key=lambda x: -barcode_reads[x])

    if not quiet:
        logger.info(f"Clustering {len(sorted_barcodes):,} barcodes "
                     f"(max edit distance {max_edit_distance}, "
                     f"min Jaccard {min_jaccard}, "
                     f"length-scaled thresholds, centroid-only Jaccard)...")

    clusters = {}           # centroid -> [member barcodes]
    barcode_to_cluster = {} # barcode -> centroid
    neighbor_to_centroid = {}  # neighbor_seq -> centroid
    
    # Centroid variant sets: fixed at creation, never updated by member merges.
    # This prevents variant set drift where successive merges gradually shift
    # the cluster's identity away from the centroid's true variant profile.
    centroid_variant_sets = {}  # centroid -> frozenset of variant_ids

    def _get_variant_set(bc):
        """Get the set of variant IDs across all reads for a barcode."""
        var_ids = set()
        if bc in barcode_read_variants:
            for read_vars in barcode_read_variants[bc]:
                for v in read_vars:
                    var_ids.add(v)
        return var_ids

    def _variants_compatible(bc_vars, centroid_vars):
        """Check if a barcode's variants are compatible with a centroid's variants."""
        bc_is_wt = len(bc_vars) == 0
        centroid_is_wt = len(centroid_vars) == 0
        # WT only matches WT
        if bc_is_wt and centroid_is_wt:
            return True
        if bc_is_wt != centroid_is_wt:
            return False
        # Both non-empty: compute Jaccard
        intersection = len(bc_vars & centroid_vars)
        union = len(bc_vars | centroid_vars)
        if union == 0:
            return True
        return (intersection / union) >= min_jaccard

    def _effective_max_dist(bc_len):
        """Scale edit distance threshold by barcode length."""
        if bc_len < expected_barcode_length:
            # Short barcodes already have errors; be more conservative
            return min(1, max_edit_distance)
        return max_edit_distance

    for i, bc in enumerate(sorted_barcodes):
        if not quiet and (i + 1) % 10000 == 0:
            logger.info(f"  Processed {i+1:,}/{len(sorted_barcodes):,} barcodes, {len(clusters):,} clusters...")

        eff_max_dist = _effective_max_dist(len(bc))
        matched_centroid = None
        bc_vars = _get_variant_set(bc)

        if eff_max_dist >= 1:
            # Check distance 0 and 1 via hash lookup
            if bc in neighbor_to_centroid:
                candidate = neighbor_to_centroid[bc]
                if _variants_compatible(bc_vars, centroid_variant_sets.get(candidate, frozenset())):
                    matched_centroid = candidate

            # Check distance 2 via d1 neighborhood overlap
            if matched_centroid is None and eff_max_dist >= 2:
                for nb in generate_levenshtein_neighbors_d1(bc):
                    if nb in neighbor_to_centroid:
                        candidate = neighbor_to_centroid[nb]
                        if _variants_compatible(bc_vars, centroid_variant_sets.get(candidate, frozenset())):
                            matched_centroid = candidate
                            break
        else:
            # Distance 0 only
            if bc in neighbor_to_centroid:
                candidate = neighbor_to_centroid[bc]
                if _variants_compatible(bc_vars, centroid_variant_sets.get(candidate, frozenset())):
                    matched_centroid = candidate

        if matched_centroid is not None:
            clusters[matched_centroid].append(bc)
            barcode_to_cluster[bc] = matched_centroid
            # No variant set update — centroid's set is immutable
        else:
            # New centroid: freeze its variant set at creation time
            clusters[bc] = [bc]
            barcode_to_cluster[bc] = bc
            centroid_variant_sets[bc] = frozenset(_get_variant_set(bc))
            if bc not in neighbor_to_centroid:
                neighbor_to_centroid[bc] = bc
            if max_edit_distance >= 1:
                for nb in generate_levenshtein_neighbors_d1(bc):
                    if nb not in neighbor_to_centroid:
                        neighbor_to_centroid[nb] = bc

    if not quiet:
        n_merged = len(sorted_barcodes) - len(clusters)
        logger.info(f"Clustering complete: {len(sorted_barcodes):,} -> {len(clusters):,} clusters ({n_merged:,} merged)")

    return clusters, barcode_to_cluster


# ---------------------------------------------------------------------------
# Position-level consensus variant calling
# ---------------------------------------------------------------------------

def consensus_variant_calling_for_cluster(cluster_barcodes, barcode_read_variants,
                                           variant_frequency_threshold=0.6,
                                           min_indel_reads=2):
    """Build position-level consensus variant call for a barcode cluster."""
    all_read_variants = []
    total_reads = 0
    for bc in cluster_barcodes:
        if bc in barcode_read_variants:
            for read_vars in barcode_read_variants[bc]:
                all_read_variants.append(read_vars)
                total_reads += 1

    if total_reads == 0:
        return "WT", 0.0, True, {}

    variant_read_count = Counter()
    for read_vars in all_read_variants:
        for var_id in read_vars:
            variant_read_count[var_id] += 1

    consensus_variants = []
    position_details = {}
    for var_id, count in variant_read_count.items():
        freq = count / total_reads
        is_indel = ('+' in var_id) or (':' in var_id and '>' not in var_id)
        passes_frequency = freq >= variant_frequency_threshold
        passes_absolute = count >= min_indel_reads if is_indel else True
        position_details[var_id] = {
            'count': count, 'total_reads': total_reads, 'frequency': freq,
            'is_indel': is_indel, 'passes_frequency': passes_frequency,
            'passes_absolute': passes_absolute, 'retained': passes_frequency and passes_absolute
        }
        if passes_frequency and passes_absolute:
            consensus_variants.append(var_id)

    consensus_tag = "WT" if not consensus_variants else json.dumps(sorted(consensus_variants))
    consensus_set = set(consensus_variants)
    matching_reads = sum(1 for rv in all_read_variants if set(rv) == consensus_set)
    support_fraction = matching_reads / total_reads if total_reads > 0 else 0.0
    is_ambiguous = support_fraction < 0.5 and total_reads > 1

    return consensus_tag, support_fraction, is_ambiguous, position_details


def resolve_barcode_clusters(clusters, barcode_variant_map, barcode_read_variants,
                              variant_frequency_threshold=0.6,
                              min_indel_reads=2, quiet=False):
    """Resolve variant assignment for each barcode cluster."""
    if not quiet:
        logger.info(f"Resolving variant consensus for {len(clusters):,} clusters...")

    cluster_resolutions = {}
    cstats = ClusteringStats(
        total_barcodes_before=sum(len(m) for m in clusters.values()),
        total_clusters_after=len(clusters)
    )
    cstats.barcodes_merged = cstats.total_barcodes_before - cstats.total_clusters_after

    for centroid, members in clusters.items():
        total_reads = sum(sum(barcode_variant_map[bc].values()) for bc in members if bc in barcode_variant_map)
        consensus_tag, support, is_ambiguous, pos_details = \
            consensus_variant_calling_for_cluster(members, barcode_read_variants, variant_frequency_threshold, min_indel_reads)
        is_low_conf = len(members) == 1 and total_reads <= 1
        cluster_resolutions[centroid] = BarcodeClusterResolution(
            consensus_variant_tag=consensus_tag, cluster_barcode=centroid,
            cluster_size=len(members), total_reads=total_reads,
            is_ambiguous=is_ambiguous, variant_support=support,
            is_low_confidence=is_low_conf, member_barcodes=members,
            position_details=pos_details)
        cstats.cluster_size_distribution[len(members)] += 1
        if is_low_conf:
            cstats.low_confidence_clusters += 1
        if is_ambiguous:
            cstats.ambiguous_clusters += 1
        else:
            cstats.resolved_clusters += 1

    cstats.singleton_clusters = cstats.cluster_size_distribution.get(1, 0)
    if not quiet:
        logger.info(f"Cluster resolution: {cstats.resolved_clusters:,} resolved, "
                     f"{cstats.ambiguous_clusters:,} ambiguous, "
                     f"{cstats.low_confidence_clusters:,} low-confidence singletons")
    return cluster_resolutions, cstats


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def variants_to_tag(variants):
    if not variants:
        return "WT"
    return json.dumps(sorted([v.id for v in variants]))


def variant_ids_from_variants(variants):
    return [v.id for v in variants]


def calculate_expected_coverage(stats, promoter_length, ref_seq, variant_types=None):
    if variant_types is None:
        variant_types = ['snv', 'ins', 'del']
    expected_snvs = promoter_length * 3 if 'snv' in variant_types else 0
    expected_insertions = promoter_length * 4 if 'ins' in variant_types else 0
    expected_deletions = promoter_length if 'del' in variant_types else 0
    observed_snvs = len(stats.observed_snvs) if 'snv' in variant_types else 0
    observed_insertions = sum(1 for v in stats.observed_insertions if '+' in v and len(v.split('+')[1]) == 1) if 'ins' in variant_types else 0
    observed_deletions = sum(1 for v in stats.observed_deletions if ':' in v and v.split(':')[1][0].isdigit() and int(''.join(c for c in v.split(':')[1] if c.isdigit())) == 1) if 'del' in variant_types else 0
    return ExpectedCoverage(
        promoter_length=promoter_length, expected_snvs=expected_snvs, observed_snvs=observed_snvs,
        snv_coverage=(observed_snvs / expected_snvs * 100) if expected_snvs > 0 else 0,
        expected_insertions=expected_insertions, observed_insertions=observed_insertions,
        insertion_coverage=(observed_insertions / expected_insertions * 100) if expected_insertions > 0 else 0,
        expected_deletions=expected_deletions, observed_deletions=observed_deletions,
        deletion_coverage=(observed_deletions / expected_deletions * 100) if expected_deletions > 0 else 0)


# ---------------------------------------------------------------------------
# BAM processing (two-pass)
# ---------------------------------------------------------------------------

def process_bam_pass1(input_bam, reference_fasta, promoter_chrom, promoter_start, promoter_end,
                      barcode_start=None, barcode_end=None, min_base_quality=20,
                      min_barcode_length=13, max_barcode_length=17,
                      extract_barcode=True, snv_only=False, quiet=False):
    """First pass: read BAM, call variants, extract barcodes."""
    stats = ProcessingStats()
    promoter_start_0 = promoter_start - 1
    promoter_end_0 = promoter_end - 1

    if extract_barcode:
        if barcode_start is None or barcode_end is None:
            raise ValueError("barcode_start and barcode_end required when extract_barcode=True")
        barcode_start_0 = barcode_start - 1
        barcode_end_0 = barcode_end
        expected_barcode_length = barcode_end - barcode_start + 1
    else:
        barcode_start_0 = None
        barcode_end_0 = None
        expected_barcode_length = 15

    bam_in = pysam.AlignmentFile(input_bam, "rb")
    fasta = pysam.FastaFile(reference_fasta)

    if promoter_chrom not in fasta.references:
        raise ValueError(f"Chromosome '{promoter_chrom}' not found in reference FASTA")
    chrom_length = fasta.get_reference_length(promoter_chrom)
    if promoter_start <= 0 or promoter_end <= 0:
        raise ValueError("Promoter coordinates must be positive (1-based)")
    if extract_barcode and (barcode_start <= 0 or barcode_end <= 0):
        raise ValueError("Barcode coordinates must be positive (1-based)")
    if promoter_end > chrom_length:
        raise ValueError(f"Promoter coordinates exceed chromosome length ({chrom_length})")
    if extract_barcode and barcode_end > chrom_length:
        raise ValueError(f"Barcode coordinates exceed chromosome length ({chrom_length})")

    ref_seq = fasta.fetch(promoter_chrom, promoter_start_0, promoter_end_0 + 1)
    if not quiet:
        logger.info(f"Pass 1: Reading BAM: {input_bam}")
        logger.info(f"Promoter: {promoter_chrom}:{promoter_start}-{promoter_end} ({promoter_end - promoter_start + 1} bp)")
        if extract_barcode:
            logger.info(f"Barcode: {promoter_chrom}:{barcode_start}-{barcode_end} ({expected_barcode_length} bp)")
        if snv_only:
            logger.info("SNV-only mode: excluding reads with indel variants in promoter")

    read_data = {}
    barcode_read_variants = defaultdict(list)

    if extract_barcode:
        fetch_start = min(promoter_start_0, barcode_start_0)
        fetch_end = max(promoter_end_0 + 1, barcode_end_0)
    else:
        fetch_start = promoter_start_0
        fetch_end = promoter_end_0 + 1

    try:
        for read in bam_in.fetch(promoter_chrom, fetch_start, fetch_end):
            stats.total_reads += 1
            if read.is_unmapped:
                stats.exclusion_counts[ExclusionReason.UNMAPPED] += 1
                continue
            if read.is_secondary:
                stats.exclusion_counts[ExclusionReason.SECONDARY] += 1
                continue
            if read.is_supplementary:
                stats.exclusion_counts[ExclusionReason.SUPPLEMENTARY] += 1
                continue

            # Single CIGAR walk for both variant calling and barcode extraction
            variants, var_status, barcode, barcode_qual, bc_status = \
                parse_promoter_and_barcode(
                    read, ref_seq, promoter_start_0, promoter_end_0,
                    bc_start_0=barcode_start_0 if extract_barcode else None,
                    bc_end_0=barcode_end_0 if extract_barcode else None,
                    min_quality=min_base_quality,
                    min_bc_length=min_barcode_length,
                    max_bc_length=max_barcode_length)
            if var_status != ExclusionReason.PASS:
                stats.exclusion_counts[var_status] += 1
                continue

            # SNV-only filter: exclude reads with any indel variant
            if snv_only and variants:
                has_indel = any(v.var_type in (VariantType.INSERTION, VariantType.DELETION) for v in variants)
                if has_indel:
                    stats.exclusion_counts[ExclusionReason.HAS_INDEL_VARIANT] += 1
                    continue

            if extract_barcode and bc_status != ExclusionReason.PASS:
                stats.exclusion_counts[bc_status] += 1
                continue

            nucleosome_count = read.get_tag("nc") if read.has_tag("nc") else None
            variant_tag = variants_to_tag(variants)
            var_ids = variant_ids_from_variants(variants)
            stats.passed_reads += 1
            num_variants = len(variants)
            stats.variants_per_read_hist[num_variants] += 1

            if num_variants == 0:
                stats.wt_reads += 1
            elif num_variants == 1:
                stats.single_variant_reads += 1
                if variants[0].var_type == VariantType.SNV:
                    stats.single_variant_snvs.add(variants[0].id)
            else:
                stats.multi_variant_reads += 1

            for var in variants:
                stats.variant_counts[var.id] += 1
                stats.variant_type_counts[var.var_type.value] += 1
                if var.var_type == VariantType.SNV:
                    stats.observed_snvs.add(var.id)
                    stats.snv_position_counts[var.position] += 1
                    if var.is_transition:
                        stats.transition_count += 1
                    else:
                        stats.transversion_count += 1
                elif var.var_type == VariantType.INSERTION:
                    stats.observed_insertions.add(var.id)
                    stats.ins_position_counts[var.position] += 1
                elif var.var_type == VariantType.DELETION:
                    stats.observed_deletions.add(var.id)
                    stats.del_position_counts[var.position] += 1
            if extract_barcode and barcode:
                stats.barcode_lengths[len(barcode)] += 1
                stats.barcode_variant_map[barcode][variant_tag] += 1
                barcode_read_variants[barcode].append(var_ids)

            if nucleosome_count is not None:
                stats.variant_nucleosome_counts[variant_tag][nucleosome_count] += 1

            # Store as tuple instead of dict (~40% less memory per read)
            # Fields: (variant_tag, var_ids, num_variants, barcode,
            #          barcode_qual, nucleosome_count, min_var_qual)
            min_var_qual = min((v.quality for v in variants if v.quality is not None), default=None)
            read_data[read.query_name] = (
                variant_tag, var_ids, num_variants, barcode,
                barcode_qual, nucleosome_count, min_var_qual,
            )

            if not quiet and stats.total_reads % 100000 == 0:
                logger.info(f"  Processed {stats.total_reads:,} reads, {stats.passed_reads:,} passed...")
    except ValueError as e:
        logger.warning(f"Could not fetch region: {e}")

    bam_in.close()
    fasta.close()
    if not quiet:
        logger.info(f"Pass 1 complete: {stats.total_reads:,} reads, {stats.passed_reads:,} passed")
        if snv_only:
            logger.info(f"  Excluded {stats.exclusion_counts.get(ExclusionReason.HAS_INDEL_VARIANT, 0):,} reads with indel variants")
    return stats, ref_seq, read_data, barcode_read_variants


def process_bam_pass2(input_bam, output_bam, promoter_chrom, promoter_start, promoter_end,
                      barcode_start, barcode_end, read_data, barcode_to_cluster,
                      cluster_resolutions, extract_barcode=True, quiet=False):
    """Second pass: write BAM with consensus-corrected tags.
    
    Returns (written, corrected, cluster_corrections) where cluster_corrections
    is a dict of centroid -> ClusterCorrection with per-cluster correction details.
    """
    if not quiet:
        logger.info(f"Pass 2: Writing consensus-corrected tags to {output_bam}")

    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    promoter_start_0 = promoter_start - 1
    promoter_end_0 = promoter_end - 1
    if extract_barcode:
        fetch_start = min(promoter_start_0, barcode_start - 1)
        fetch_end = max(promoter_end_0 + 1, barcode_end)
    else:
        fetch_start = promoter_start_0
        fetch_end = promoter_end_0 + 1

    written = 0
    corrected = 0

    # Per-cluster correction tracking:
    # cluster_raw_tags[centroid] = Counter of raw_tag -> count
    # cluster_read_counts[centroid] = total reads in cluster
    cluster_raw_tags = defaultdict(Counter)
    cluster_read_counts = defaultdict(int)
    cluster_corrected_counts = defaultdict(int)
    # Per-cluster, per-variant tracking: how many reads carried each variant
    cluster_raw_variant_counts = defaultdict(Counter)  # centroid -> variant_id -> read_count

    # read_data tuple indices:
    #   0=variant_tag, 1=var_ids, 2=num_variants, 3=barcode,
    #   4=barcode_qual, 5=nucleosome_count, 6=min_var_qual
    RD_VARIANT_TAG = 0
    RD_VAR_IDS = 1
    RD_BARCODE = 3
    RD_BARCODE_QUAL = 4
    RD_MIN_VAR_QUAL = 6

    try:
        for read in bam_in.fetch(promoter_chrom, fetch_start, fetch_end):
            # Skip supplementary alignments — they share query_name with
            # their primary and would incorrectly inherit its tags
            if read.is_supplementary or read.is_secondary:
                continue
            if read.query_name not in read_data:
                continue
            rd = read_data[read.query_name]
            raw_tag = rd[RD_VARIANT_TAG]
            barcode = rd[RD_BARCODE]

            if extract_barcode and barcode and barcode in barcode_to_cluster:
                centroid = barcode_to_cluster[barcode]
                if centroid in cluster_resolutions:
                    res = cluster_resolutions[centroid]
                    consensus_tag = res.consensus_variant_tag
                    cluster_size = res.cluster_size
                else:
                    consensus_tag = raw_tag
                    cluster_size = 1
            else:
                consensus_tag = raw_tag
                centroid = barcode
                cluster_size = 1

            if consensus_tag == "WT":
                n_cons_var = 0
            else:
                try:
                    n_cons_var = len(json.loads(consensus_tag))
                except (json.JSONDecodeError, TypeError):
                    n_cons_var = 0

            read.set_tag("PV", consensus_tag)
            read.set_tag("VC", n_cons_var)
            read.set_tag("PR", raw_tag)
            if extract_barcode and barcode:
                read.set_tag("BC", barcode)
                if rd[RD_BARCODE_QUAL] is not None:
                    read.set_tag("BQ", int(round(rd[RD_BARCODE_QUAL])))
                if centroid:
                    read.set_tag("BK", centroid)
                read.set_tag("CS", cluster_size)
            if rd[RD_MIN_VAR_QUAL] is not None:
                read.set_tag("VQ", rd[RD_MIN_VAR_QUAL])

            bam_out.write(read)
            written += 1

            # Track correction details per cluster
            if centroid is not None:
                cluster_raw_tags[centroid][raw_tag] += 1
                cluster_read_counts[centroid] += 1
                if consensus_tag != raw_tag:
                    corrected += 1
                    cluster_corrected_counts[centroid] += 1
                # Track per-variant read counts for this cluster
                for var_id in rd[RD_VAR_IDS]:
                    cluster_raw_variant_counts[centroid][var_id] += 1

    except ValueError as e:
        logger.warning(f"Could not fetch region in pass 2: {e}")

    bam_in.close()
    bam_out.close()

    # Build per-cluster correction summaries (only for clusters with corrections)
    cluster_corrections = {}
    for centroid, corr_count in cluster_corrected_counts.items():
        if corr_count == 0:
            continue
        total_reads = cluster_read_counts[centroid]
        raw_tag_counts = cluster_raw_tags[centroid]
        most_common_raw = raw_tag_counts.most_common(1)[0][0]

        # Get consensus variant set
        consensus_tag = "WT"
        if centroid in cluster_resolutions:
            consensus_tag = cluster_resolutions[centroid].consensus_variant_tag
        consensus_vars = set()
        if consensus_tag != "WT":
            try:
                consensus_vars = set(json.loads(consensus_tag))
            except (json.JSONDecodeError, TypeError):
                pass

        # Get all raw variant IDs seen in this cluster
        raw_var_counts = cluster_raw_variant_counts[centroid]
        all_raw_vars = set(raw_var_counts.keys())

        # Variants removed: appeared in raw calls but not in consensus
        variants_removed = sorted(all_raw_vars - consensus_vars)
        # Variants added: in consensus but some reads didn't have them
        variants_added = sorted(consensus_vars - all_raw_vars) if consensus_vars else []

        # Frequency of each removed variant (fraction of cluster reads that had it)
        removed_freqs = {}
        for var_id in variants_removed:
            removed_freqs[var_id] = raw_var_counts[var_id] / total_reads if total_reads > 0 else 0.0

        cluster_size = cluster_resolutions[centroid].cluster_size if centroid in cluster_resolutions else 1
        cluster_corrections[centroid] = ClusterCorrection(
            centroid=centroid, cluster_size=cluster_size,
            total_reads=total_reads, raw_variant_tag=most_common_raw,
            consensus_variant_tag=consensus_tag,
            reads_corrected=corr_count,
            correction_rate=corr_count / total_reads if total_reads > 0 else 0.0,
            variants_removed=variants_removed, variants_added=variants_added,
            removed_variant_freqs=removed_freqs)

    if not quiet:
        logger.info(f"Pass 2 complete: wrote {written:,} reads, {corrected:,} corrected by consensus")
        logger.info(f"  {len(cluster_corrections):,} clusters had corrections applied")
    return written, corrected, cluster_corrections


# ---------------------------------------------------------------------------
# Report writing
# ---------------------------------------------------------------------------

def write_stats_report(stats, output_path, clustering_stats=None):
    with open(output_path, 'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"total_reads\t{stats.total_reads}\n")
        f.write(f"passed_reads\t{stats.passed_reads}\n")
        f.write(f"pass_rate\t{stats.passed_reads / stats.total_reads:.4f}\n" if stats.total_reads > 0 else "pass_rate\t0\n")
        for reason in ExclusionReason:
            if reason != ExclusionReason.PASS:
                f.write(f"excluded_{reason.value}\t{stats.exclusion_counts.get(reason, 0)}\n")
        f.write(f"wt_reads\t{stats.wt_reads}\n")
        f.write(f"single_variant_reads\t{stats.single_variant_reads}\n")
        f.write(f"multi_variant_reads\t{stats.multi_variant_reads}\n")
        f.write(f"unique_variants\t{len(stats.variant_counts)}\n")
        f.write(f"unique_barcodes\t{len(stats.barcode_variant_map)}\n")
        ti_tv = stats.transition_count / stats.transversion_count if stats.transversion_count > 0 else float('inf')
        f.write(f"ti_tv_ratio\t{ti_tv:.4f}\n")
        if clustering_stats:
            f.write(f"barcode_clusters\t{clustering_stats.total_clusters_after}\n")
            f.write(f"barcodes_merged\t{clustering_stats.barcodes_merged}\n")
            f.write(f"singleton_clusters\t{clustering_stats.singleton_clusters}\n")
            f.write(f"resolved_clusters\t{clustering_stats.resolved_clusters}\n")
            f.write(f"ambiguous_clusters\t{clustering_stats.ambiguous_clusters}\n")
            f.write(f"low_confidence_clusters\t{clustering_stats.low_confidence_clusters}\n")


def write_variants_report(stats, output_path):
    variant_to_barcodes = defaultdict(set)
    for barcode, var_counts in stats.barcode_variant_map.items():
        for var_tag in var_counts.keys():
            variant_to_barcodes[var_tag].add(barcode)
    with open(output_path, 'w') as f:
        f.write("variant_id\tread_count\tbarcode_count\tvariant_type\tis_transition\n")
        for var_id, count in sorted(stats.variant_counts.items(), key=lambda x: -x[1]):
            if ">" in var_id:
                var_type = "snv"
                parts = var_id.split(":")
                ref, alt = parts[1].split(">") if len(parts) == 2 and ">" in parts[1] else ("", "")
                is_trans = (ref.upper(), alt.upper()) in {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
            elif "+" in var_id:
                var_type, is_trans = "insertion", False
            else:
                var_type, is_trans = "deletion", False
            single_var_tag = json.dumps([var_id])
            bc_count = len(variant_to_barcodes.get(single_var_tag, set()))
            f.write(f"{var_id}\t{count}\t{bc_count}\t{var_type}\t{is_trans}\n")


def write_cluster_resolution_report(cluster_resolutions, output_path):
    with open(output_path, 'w') as f:
        f.write("cluster_barcode\tconsensus_variant\tcluster_size\ttotal_reads\t"
                "variant_support\tis_ambiguous\tis_low_confidence\tmember_barcodes\n")
        for centroid in sorted(cluster_resolutions.keys()):
            res = cluster_resolutions[centroid]
            members_str = ",".join(res.member_barcodes) if len(res.member_barcodes) <= 10 else \
                          ",".join(res.member_barcodes[:10]) + f"...+{len(res.member_barcodes)-10}more"
            f.write(f"{centroid}\t{res.consensus_variant_tag}\t{res.cluster_size}\t"
                    f"{res.total_reads}\t{res.variant_support:.4f}\t{res.is_ambiguous}\t"
                    f"{res.is_low_confidence}\t{members_str}\n")


def write_exclusion_report(stats, output_path):
    with open(output_path, 'w') as f:
        f.write("exclusion_reason\tcount\tproportion\n")
        for reason in ExclusionReason:
            count = stats.passed_reads if reason == ExclusionReason.PASS else stats.exclusion_counts.get(reason, 0)
            prop = count / stats.total_reads if stats.total_reads > 0 else 0
            f.write(f"{reason.value}\t{count}\t{prop:.6f}\n")


def write_expected_coverage_report(coverage, output_path):
    with open(output_path, 'w') as f:
        f.write("variant_type\texpected\tobserved\tcoverage_percent\n")
        f.write(f"snv\t{coverage.expected_snvs}\t{coverage.observed_snvs}\t{coverage.snv_coverage:.2f}\n")
        f.write(f"single_bp_insertion\t{coverage.expected_insertions}\t{coverage.observed_insertions}\t{coverage.insertion_coverage:.2f}\n")
        f.write(f"single_bp_deletion\t{coverage.expected_deletions}\t{coverage.observed_deletions}\t{coverage.deletion_coverage:.2f}\n")


def write_correction_report(cluster_corrections, output_path):
    """Write per-cluster correction details as TSV."""
    with open(output_path, 'w') as f:
        f.write("cluster_centroid\tcluster_size\ttotal_reads\tmost_common_raw_tag\t"
                "consensus_tag\treads_corrected\tcorrection_rate\t"
                "variants_removed\tremoved_variant_frequencies\tvariants_added\n")
        for centroid in sorted(cluster_corrections.keys()):
            cc = cluster_corrections[centroid]
            removed_str = ",".join(cc.variants_removed) if cc.variants_removed else "none"
            removed_freq_str = ",".join(
                f"{v}:{cc.removed_variant_freqs[v]:.3f}" for v in cc.variants_removed
            ) if cc.variants_removed else "none"
            added_str = ",".join(cc.variants_added) if cc.variants_added else "none"
            f.write(f"{centroid}\t{cc.cluster_size}\t{cc.total_reads}\t"
                    f"{cc.raw_variant_tag}\t{cc.consensus_variant_tag}\t"
                    f"{cc.reads_corrected}\t{cc.correction_rate:.4f}\t"
                    f"{removed_str}\t{removed_freq_str}\t{added_str}\n")


# ---------------------------------------------------------------------------
# Plot generation
# ---------------------------------------------------------------------------

def generate_plots(stats, cluster_resolutions=None, clustering_stats=None,
                   coverage=None, cluster_corrections=None,
                   promoter_start=0, promoter_end=0,
                   expected_barcode_length=15,
                   top_variants_nucleosome=10, output_prefix=None, save_individual_pngs=False):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available, skipping plot generation")
        return {}

    figures = {}
    plt.rcParams['figure.figsize'] = (10, 6)
    plt.rcParams['font.size'] = 12

    def _save(fig, name):
        if save_individual_pngs and output_prefix:
            fig.savefig(f"{output_prefix}_{name}.png", dpi=150, bbox_inches='tight')

    # 1. Variants per read (from pre-binned Counter)
    if stats.variants_per_read_hist:
        fig, ax = plt.subplots()
        max_var = max(stats.variants_per_read_hist.keys())
        x_vals = list(range(0, min(max_var + 2, 20)))
        y_vals = [stats.variants_per_read_hist.get(x, 0) for x in x_vals]
        ax.bar(x_vals, y_vals, edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xlabel('Number of Variants per Read')
        ax.set_ylabel('Read Count')
        ax.set_title('Distribution of Variants per Read')
        total = sum(stats.variants_per_read_hist.values())
        ax.text(0.95, 0.95, f'WT: {stats.wt_reads/total*100:.1f}%\nSingle: {stats.single_variant_reads/total*100:.1f}%',
                transform=ax.transAxes, ha='right', va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        figures['variants_per_read'] = fig
        _save(fig, 'variants_per_read')

    # 2. Barcode length distribution
    if stats.barcode_lengths:
        fig, ax = plt.subplots()
        lengths = sorted(stats.barcode_lengths.keys())
        counts = [stats.barcode_lengths[l] for l in lengths]
        colors = ['forestgreen' if l == expected_barcode_length else 'steelblue' for l in lengths]
        ax.bar(lengths, counts, color=colors, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Barcode Length (bp)')
        ax.set_ylabel('Read Count')
        ax.set_title('Barcode Length Distribution')
        total_bc = sum(counts)
        expected_count = stats.barcode_lengths.get(expected_barcode_length, 0)
        ax.axvline(expected_barcode_length, color='red', linestyle='--', alpha=0.5)
        ax.text(0.95, 0.95,
                f'Expected ({expected_barcode_length}bp): {expected_count:,} ({expected_count/total_bc*100:.1f}%)\n'
                f'Total: {total_bc:,}',
                transform=ax.transAxes, ha='right', va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        figures['barcode_length_distribution'] = fig
        _save(fig, 'barcode_length_distribution')

    # 3. Cluster resolution status + cluster size distribution
    if cluster_resolutions and clustering_stats:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        ax = axes[0]
        resolved = clustering_stats.resolved_clusters
        ambiguous = clustering_stats.ambiguous_clusters
        low_conf = clustering_stats.low_confidence_clusters
        ax.bar(['Resolved', 'Ambiguous'], [resolved, ambiguous], color=['forestgreen', 'firebrick'], alpha=0.7, edgecolor='black')
        for i, v in enumerate([resolved, ambiguous]):
            ax.text(i, v + max(resolved, ambiguous) * 0.02, f'{v:,}', ha='center', va='bottom', fontweight='bold')
        ax.set_ylabel('Number of Barcode Clusters')
        ax.set_title('Cluster Resolution\n(Edit-Distance Clustering + Position-Level Consensus)')
        ax.text(0.95, 0.95, f'Low-confidence\nsingletons: {low_conf:,}',
                transform=ax.transAxes, ha='right', va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

        ax = axes[1]
        size_dist = clustering_stats.cluster_size_distribution
        sizes = sorted(size_dist.keys())
        if max(sizes) > 20:
            display_sizes = [s for s in sizes if s <= 19]
            display_counts = [size_dist[s] for s in display_sizes]
            display_sizes.append(20)
            display_counts.append(sum(size_dist[s] for s in sizes if s >= 20))
            xlabels = [str(s) for s in display_sizes[:-1]] + ['20+']
        else:
            display_sizes = sizes
            display_counts = [size_dist[s] for s in sizes]
            xlabels = [str(s) for s in display_sizes]
        ax.bar(range(len(display_sizes)), display_counts, color='steelblue', alpha=0.7, edgecolor='black')
        ax.set_xticks(range(len(display_sizes)))
        ax.set_xticklabels(xlabels)
        ax.set_xlabel('Barcodes per Cluster')
        ax.set_ylabel('Number of Clusters')
        ax.set_title('Cluster Size Distribution')
        ax.set_yscale('log')
        plt.tight_layout()
        figures['cluster_resolution'] = fig
        _save(fig, 'cluster_resolution')

    # 4. Post-clustering: barcode clusters per consensus variant (top 50 + WT)
    if cluster_resolutions:
        variant_to_clusters = defaultdict(int)
        variant_to_reads = defaultdict(int)
        for res in cluster_resolutions.values():
            if not res.is_ambiguous:
                variant_to_clusters[res.consensus_variant_tag] += 1
                variant_to_reads[res.consensus_variant_tag] += res.total_reads

        if variant_to_clusters:
            # Sort by cluster count, put WT first if present
            sorted_variants = sorted(variant_to_clusters.items(), key=lambda x: -x[1])
            # Separate WT
            wt_entry = None
            non_wt = []
            for tag, count in sorted_variants:
                if tag == "WT":
                    wt_entry = (tag, count)
                else:
                    non_wt.append((tag, count))

            # Take top 50 non-WT
            display_items = non_wt[:50]
            if wt_entry:
                display_items = [wt_entry] + display_items

            fig, ax = plt.subplots(figsize=(14, 6))
            labels = []
            cluster_counts = []
            read_counts = []
            for tag, cc in display_items:
                if tag == "WT":
                    labels.append("WT")
                else:
                    try:
                        var_list = json.loads(tag)
                        labels.append(", ".join(var_list) if len(var_list) <= 2 else f"{var_list[0]}... +{len(var_list)-1}")
                    except (json.JSONDecodeError, TypeError):
                        labels.append(tag[:20])
                cluster_counts.append(cc)
                read_counts.append(variant_to_reads.get(tag, 0))

            x = np.arange(len(labels))
            ax.bar(x, cluster_counts, color='steelblue', edgecolor='black', alpha=0.7)
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=90, fontsize=7)
            ax.set_xlabel('Consensus Variant')
            ax.set_ylabel('Number of Barcode Clusters')
            ax.set_title('Barcode Clusters per Variant (Post-Clustering, Resolved Only)')
            # Add read count as secondary annotation for top entries
            for i in range(min(5, len(x))):
                ax.text(i, cluster_counts[i] + max(cluster_counts) * 0.01,
                       f'{read_counts[i]:,}r', ha='center', va='bottom', fontsize=7, color='gray')
            ax.text(0.95, 0.95, f'Total resolved: {sum(cluster_counts):,} clusters\n'
                    f'Unique variants: {len(variant_to_clusters):,}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            plt.tight_layout()
            figures['clusters_per_variant'] = fig
            _save(fig, 'clusters_per_variant')

    # 5-7. Variant read counts by type
    for var_filter, color, label in [
        (lambda k: '>' in k, 'steelblue', 'SNV'),
        (lambda k: '+' in k, 'darkorange', 'Insertion'),
        (lambda k: '>' not in k and '+' not in k, 'firebrick', 'Deletion')
    ]:
        filtered = {k: v for k, v in stats.variant_counts.items() if var_filter(k)}
        if filtered:
            fig, ax = plt.subplots(figsize=(12, 6))
            sorted_vars = sorted(filtered.items(), key=lambda x: -x[1])[:50]
            ax.bar(range(len(sorted_vars)), [v[1] for v in sorted_vars], edgecolor='black', alpha=0.7, color=color)
            ax.set_xlabel(f'{label} Variant')
            ax.set_ylabel('Read Count')
            ax.set_title(f'Reads per {label} Variant (Top {len(sorted_vars)})')
            ax.set_xticks(range(len(sorted_vars)))
            ax.set_xticklabels([v[0] for v in sorted_vars], rotation=90, fontsize=8)
            plt.tight_layout()
            figures[f'{label.lower()}_read_counts'] = fig
            _save(fig, f'{label.lower()}_read_counts')

    # 8. Position coverage
    if stats.snv_position_counts or stats.ins_position_counts or stats.del_position_counts:
        fig, ax = plt.subplots(figsize=(14, 6))
        positions = list(range(promoter_start, promoter_end + 1))
        snv_counts = [stats.snv_position_counts.get(p, 0) for p in positions]
        ins_counts = [stats.ins_position_counts.get(p, 0) for p in positions]
        del_counts = [stats.del_position_counts.get(p, 0) for p in positions]
        ax.bar(positions, snv_counts, label='SNVs', alpha=0.7, color='steelblue')
        ax.bar(positions, ins_counts, bottom=snv_counts, label='Insertions', alpha=0.7, color='darkorange')
        ax.bar(positions, del_counts, bottom=[s+i for s, i in zip(snv_counts, ins_counts)], label='Deletions', alpha=0.7, color='firebrick')
        ax.set_xlabel('Position in Promoter')
        ax.set_ylabel('Variant Count')
        ax.set_title('Variant Frequency by Position')
        ax.legend()
        plt.tight_layout()
        figures['position_coverage'] = fig
        _save(fig, 'position_coverage')

    # 9. Position frequency
    if stats.passed_reads > 0:
        fig, ax = plt.subplots(figsize=(14, 6))
        positions = list(range(promoter_start, promoter_end + 1))
        total_variants = [stats.snv_position_counts.get(p, 0) + stats.ins_position_counts.get(p, 0) + stats.del_position_counts.get(p, 0) for p in positions]
        frequencies = [c / stats.passed_reads * 100 for c in total_variants]
        ax.bar(positions, frequencies, alpha=0.7, color='steelblue')
        ax.set_xlabel('Position in Promoter')
        ax.set_ylabel('Variant Frequency (%)')
        ax.set_title('Per-Position Variant Frequency')
        ax.axhline(np.mean(frequencies), color='red', linestyle='--', label=f'Mean: {np.mean(frequencies):.2f}%')
        ax.legend()
        plt.tight_layout()
        figures['position_frequency'] = fig
        _save(fig, 'position_frequency')

    # 10. Expected coverage
    if coverage:
        fig, ax = plt.subplots(figsize=(8, 6))
        categories = ['SNVs', 'Single-bp\nInsertions', 'Single-bp\nDeletions']
        expected = [coverage.expected_snvs, coverage.expected_insertions, coverage.expected_deletions]
        observed = [coverage.observed_snvs, coverage.observed_insertions, coverage.observed_deletions]
        x = np.arange(len(categories))
        width = 0.35
        ax.bar(x - width/2, expected, width, label='Expected', alpha=0.7, color='lightgray', edgecolor='black')
        ax.bar(x + width/2, observed, width, label='Observed', alpha=0.7, color='steelblue', edgecolor='black')
        ax.set_ylabel('Count')
        ax.set_title('Expected vs Observed Variant Coverage')
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.legend()
        for i, (exp, obs) in enumerate(zip(expected, observed)):
            if exp > 0:
                ax.annotate(f'{obs/exp*100:.1f}%', xy=(i + width/2, obs), ha='center', va='bottom', fontsize=10)
        plt.tight_layout()
        figures['expected_coverage'] = fig
        _save(fig, 'expected_coverage')

    # 11. Nucleosome distribution (from Counter-based accumulators)
    if stats.variant_nucleosome_counts:
        variant_read_counts = {tag: sum(counter.values())
                               for tag, counter in stats.variant_nucleosome_counts.items()}
        sorted_variants = sorted(variant_read_counts.items(), key=lambda x: -x[1])
        variants_to_plot = []
        if 'WT' in stats.variant_nucleosome_counts:
            variants_to_plot.append('WT')
        for var_tag, _ in sorted_variants:
            if var_tag != 'WT' and len(variants_to_plot) < top_variants_nucleosome + 1:
                variants_to_plot.append(var_tag)
        if variants_to_plot:
            nrows = (len(variants_to_plot) + 2) // 3
            ncols = min(3, len(variants_to_plot))
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4 * nrows))
            axes = [axes] if len(variants_to_plot) == 1 else (axes.flatten() if hasattr(axes, 'flatten') else [axes])
            for idx, var_tag in enumerate(variants_to_plot):
                if idx < len(axes):
                    ax = axes[idx]
                    nc_counter = stats.variant_nucleosome_counts[var_tag]
                    if nc_counter:
                        max_nc = max(nc_counter.keys())
                        x_vals = list(range(0, max_nc + 2))
                        y_vals = [nc_counter.get(x, 0) for x in x_vals]
                        ax.bar(x_vals, y_vals, edgecolor='black', alpha=0.7, color='steelblue')
                        ax.set_xlabel('Nucleosome Count')
                        ax.set_ylabel('Read Count')
                        label = var_tag if len(var_tag) < 30 else var_tag[:27] + '...'
                        n_reads = sum(nc_counter.values())
                        ax.set_title(f'{label}\n(n={n_reads})')
                        mean_nc = sum(k * v for k, v in nc_counter.items()) / n_reads
                        ax.axvline(mean_nc, color='red', linestyle='--', label=f'Mean: {mean_nc:.1f}')
                        ax.legend(fontsize=8)
            for idx in range(len(variants_to_plot), len(axes)):
                axes[idx].set_visible(False)
            plt.suptitle('Nucleosome Count Distribution by Variant', fontsize=14, y=1.02)
            plt.tight_layout()
            figures['nucleosome_distribution'] = fig
            _save(fig, 'nucleosome_distribution')

    # 12. Consensus correction: correction rate distribution across clusters
    if cluster_corrections:
        correction_rates = [cc.correction_rate for cc in cluster_corrections.values()]
        if correction_rates:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            # Panel A: Distribution of correction rates across corrected clusters
            ax = axes[0]
            ax.hist(correction_rates, bins=20, edgecolor='black', alpha=0.7, color='steelblue',
                    range=(0, 1))
            ax.set_xlabel('Correction Rate (fraction of reads corrected)')
            ax.set_ylabel('Number of Clusters')
            ax.set_title('Per-Cluster Correction Rate\n(clusters with ≥1 correction)')
            ax.axvline(np.median(correction_rates), color='red', linestyle='--',
                       label=f'Median: {np.median(correction_rates):.2f}')
            ax.legend()
            ax.text(0.95, 0.95, f'Clusters with corrections: {len(correction_rates):,}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            # Panel B: Frequency distribution of removed variants
            # Shows at what read frequency variants are being dropped
            all_removed_freqs = []
            for cc in cluster_corrections.values():
                all_removed_freqs.extend(cc.removed_variant_freqs.values())
            ax = axes[1]
            if all_removed_freqs:
                ax.hist(all_removed_freqs, bins=20, edgecolor='black', alpha=0.7,
                        color='firebrick', range=(0, 1))
                ax.set_xlabel('Read Frequency of Removed Variant')
                ax.set_ylabel('Count (variant × cluster instances)')
                ax.set_title('Frequency of Variants Removed by Consensus')
                ax.axvline(0.6, color='black', linestyle=':', alpha=0.5,
                           label='Consensus threshold (0.6)')
                ax.axvline(np.median(all_removed_freqs), color='red', linestyle='--',
                           label=f'Median: {np.median(all_removed_freqs):.2f}')
                ax.legend(fontsize=9)
                ax.text(0.95, 0.95,
                        f'Total removals: {len(all_removed_freqs):,}\n'
                        f'Below 20%: {sum(1 for f in all_removed_freqs if f < 0.2):,}\n'
                        f'20-60%: {sum(1 for f in all_removed_freqs if 0.2 <= f < 0.6):,}',
                        transform=ax.transAxes, ha='right', va='top', fontsize=9,
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            else:
                ax.text(0.5, 0.5, 'No variants removed', ha='center', va='center',
                        transform=ax.transAxes, fontsize=14, color='gray')
                ax.set_title('Frequency of Variants Removed by Consensus')

            plt.tight_layout()
            figures['consensus_corrections'] = fig
            _save(fig, 'consensus_corrections')

    return figures


# ---------------------------------------------------------------------------
# HTML / PDF report
# ---------------------------------------------------------------------------

def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def generate_html_report(stats, cluster_resolutions, clustering_stats,
                          coverage, figures, output_path, sample_name, include_barcode=True):
    pass_rate = stats.passed_reads / stats.total_reads * 100 if stats.total_reads > 0 else 0
    ti_tv = stats.transition_count / stats.transversion_count if stats.transversion_count > 0 else float('inf')
    wt_pct = stats.wt_reads / stats.passed_reads * 100 if stats.passed_reads > 0 else 0

    barcode_card = ""
    cluster_card = ""
    if include_barcode and clustering_stats:
        barcode_card = (f'<div class="card"><h3>Unique Barcodes</h3>'
                       f'<div class="value">{clustering_stats.total_barcodes_before:,}</div>'
                       f'<div class="subvalue">&rarr; {clustering_stats.total_clusters_after:,} clusters</div></div>')
        cluster_card = (f'<div class="card"><h3>Cluster Resolution</h3>'
                       f'<div class="value">{clustering_stats.resolved_clusters:,}</div>'
                       f'<div class="subvalue">{clustering_stats.ambiguous_clusters:,} ambiguous, '
                       f'{clustering_stats.low_confidence_clusters:,} low-conf</div></div>')

    html = f"""<!DOCTYPE html>
<html><head><title>MPRA Report: {sample_name}</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; background: #f5f5f5; }}
.header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px; }}
.header h1 {{ margin: 0 0 10px 0; }} .header p {{ margin: 0; opacity: 0.9; }}
.summary-cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 30px; }}
.card {{ background: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
.card h3 {{ margin: 0 0 10px 0; color: #666; font-size: 14px; text-transform: uppercase; }}
.card .value {{ font-size: 28px; font-weight: bold; color: #333; }}
.card .subvalue {{ font-size: 14px; color: #888; }}
.section {{ background: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 30px; }}
.section h2 {{ margin-top: 0; color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px; }}
.plot-container {{ text-align: center; margin: 20px 0; }}
.plot-container img {{ max-width: 100%; border-radius: 5px; }}
.plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }}
table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
th {{ background-color: #f8f9fa; font-weight: 600; }}
.good {{ color: #28a745; }} .warning {{ color: #ffc107; }} .bad {{ color: #dc3545; }}
</style></head><body>
<div class="header"><h1>MPRA Variant-Barcode QC Report</h1><p>Sample: {sample_name}</p></div>
<div class="summary-cards">
<div class="card"><h3>Total Reads</h3><div class="value">{stats.total_reads:,}</div></div>
<div class="card"><h3>Passed Reads</h3><div class="value">{stats.passed_reads:,}</div><div class="subvalue">{pass_rate:.1f}% pass rate</div></div>
<div class="card"><h3>Unique Variants</h3><div class="value">{len(stats.variant_counts):,}</div></div>
{barcode_card}{cluster_card}
<div class="card"><h3>WT Reads</h3><div class="value">{stats.wt_reads:,}</div><div class="subvalue">{wt_pct:.1f}% of passed</div></div>
<div class="card"><h3>Ti/Tv Ratio</h3><div class="value">{ti_tv:.2f}</div><div class="subvalue">{stats.transition_count:,} Ti / {stats.transversion_count:,} Tv</div></div>
</div>"""

    # Read classification
    html += f"""<div class="section"><h2>Read Classification</h2>
<table><tr><th>Category</th><th>Count</th><th>Percentage</th></tr>
<tr><td>Wild-type reads</td><td>{stats.wt_reads:,}</td><td>{stats.wt_reads/stats.passed_reads*100:.2f}%</td></tr>
<tr><td>Single variant reads</td><td>{stats.single_variant_reads:,}</td><td>{stats.single_variant_reads/stats.passed_reads*100:.2f}%</td></tr>
<tr><td>Multi-variant reads</td><td>{stats.multi_variant_reads:,}</td><td>{stats.multi_variant_reads/stats.passed_reads*100:.2f}%</td></tr>
</table><div class="plot-grid">"""
    if 'variants_per_read' in figures:
        html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures["variants_per_read"])}"></div>'
    html += '</div></div>'

    # Expected coverage
    html += f"""<div class="section"><h2>Expected Variant Coverage</h2>
<table><tr><th>Variant Type</th><th>Expected</th><th>Observed</th><th>Coverage</th></tr>
<tr><td>SNVs</td><td>{coverage.expected_snvs:,}</td><td>{coverage.observed_snvs:,}</td>
<td class="{'good' if coverage.snv_coverage > 80 else 'warning' if coverage.snv_coverage > 50 else 'bad'}">{coverage.snv_coverage:.1f}%</td></tr>
<tr><td>Single-bp Insertions</td><td>{coverage.expected_insertions:,}</td><td>{coverage.observed_insertions:,}</td>
<td class="{'good' if coverage.insertion_coverage > 80 else 'warning' if coverage.insertion_coverage > 50 else 'bad'}">{coverage.insertion_coverage:.1f}%</td></tr>
<tr><td>Single-bp Deletions</td><td>{coverage.expected_deletions:,}</td><td>{coverage.observed_deletions:,}</td>
<td class="{'good' if coverage.deletion_coverage > 80 else 'warning' if coverage.deletion_coverage > 50 else 'bad'}">{coverage.deletion_coverage:.1f}%</td></tr>
</table>"""
    if 'expected_coverage' in figures:
        html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures["expected_coverage"])}"></div>'
    html += '</div>'

    # Variant distribution
    html += '<div class="section"><h2>Variant Distribution</h2><div class="plot-grid">'
    for pn in ['position_coverage', 'position_frequency']:
        if pn in figures:
            html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures[pn])}"></div>'
    html += '</div></div>'

    # Variant read counts
    html += '<div class="section"><h2>Reads per Variant Type</h2><div class="plot-grid">'
    for pn in ['snv_read_counts', 'insertion_read_counts', 'deletion_read_counts']:
        if pn in figures:
            html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures[pn])}"></div>'
    html += '</div></div>'

    # Barcode / Cluster section
    if include_barcode and clustering_stats:
        html += '<div class="section"><h2>Barcode Clustering &amp; Resolution</h2>'
        html += f"""<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Raw unique barcodes</td><td>{clustering_stats.total_barcodes_before:,}</td></tr>
<tr><td>Clusters after merging</td><td>{clustering_stats.total_clusters_after:,}</td></tr>
<tr><td>Barcodes merged</td><td>{clustering_stats.barcodes_merged:,}</td></tr>
<tr><td>Singleton clusters</td><td>{clustering_stats.singleton_clusters:,}</td></tr>
<tr><td>Resolved clusters</td><td class="good">{clustering_stats.resolved_clusters:,}</td></tr>
<tr><td>Ambiguous clusters</td><td class="bad">{clustering_stats.ambiguous_clusters:,}</td></tr>
<tr><td>Low-confidence singletons</td><td class="warning">{clustering_stats.low_confidence_clusters:,}</td></tr>
</table><div class="plot-grid">"""
        for pn in ['cluster_resolution', 'clusters_per_variant', 'barcode_length_distribution']:
            if pn in figures:
                html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures[pn])}"></div>'
        html += '</div></div>'

    # Nucleosome
    if 'nucleosome_distribution' in figures:
        html += f'<div class="section"><h2>Nucleosome Distribution by Variant</h2><div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures["nucleosome_distribution"])}"></div></div>'

    # Consensus corrections
    if 'consensus_corrections' in figures:
        html += '<div class="section"><h2>Consensus Corrections</h2>'
        html += '<p>Shows how consensus calling corrected per-read variant tags. '
        html += 'Left: distribution of correction rates across clusters that had at least one read corrected. '
        html += 'Right: the read frequency at which removed variants appeared in their cluster — '
        html += 'variants near 0% are likely sequencing errors, while variants near the 0.6 threshold '
        html += 'may represent borderline signal.</p>'
        html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures["consensus_corrections"])}"></div>'
        html += '</div>'

    html += '</body></html>'
    with open(output_path, 'w') as f:
        f.write(html)
    logger.info(f"Generated HTML report: {output_path}")


def generate_pdf_report(figures, output_path, sample_name):
    try:
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available, skipping PDF generation")
        return
    with PdfPages(output_path) as pdf:
        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis('off')
        ax.text(0.5, 0.6, 'MPRA Variant-Barcode QC Report', ha='center', va='center', fontsize=24, fontweight='bold')
        ax.text(0.5, 0.4, f'Sample: {sample_name}', ha='center', va='center', fontsize=16)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        for name, fig in figures.items():
            pdf.savefig(fig, bbox_inches='tight')
    logger.info(f"Generated PDF report: {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_region_coords(region_str):
    """Parse 'START-END' or 'CHROM:START-END' into (start, end) integers (1-based inclusive).

    Returns (chrom_or_None, start, end).
    """
    if ':' in region_str:
        chrom, coords = region_str.split(':', 1)
        parts = coords.split('-')
    else:
        chrom = None
        parts = region_str.split('-')
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"Region must be CHROM:START-END or START-END (1-based inclusive), got: {region_str}")
    try:
        start, end = int(parts[0]), int(parts[1])
    except ValueError:
        raise argparse.ArgumentTypeError(f"Region coordinates must be integers, got: {region_str}")
    if start <= 0 or end <= 0:
        raise argparse.ArgumentTypeError(f"Coordinates must be positive (1-based), got: {region_str}")
    if start > end:
        raise argparse.ArgumentTypeError(f"Start must be <= end, got: {region_str}")
    return chrom, start, end


def parse_args():
    parser = argparse.ArgumentParser(
        description="Tag MPRA reads with promoter variants and barcodes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s input.bam reference.fasta \\
        --promoter-chrom LDLR --promoter-region 3184-3501 \\
        --barcode-region 3294-3309

    # SNV-only mode:
    %(prog)s input.bam reference.fasta \\
        --promoter-chrom LDLR --promoter-region 3184-3501 \\
        --barcode-region 3294-3309 --snv-only

    # Without barcode extraction:
    %(prog)s input.bam reference.fasta \\
        --promoter-chrom LDLR --promoter-region 3184-3501 --no-barcode

Output BAM tags:
    PV  Consensus variant (post-clustering)
    PR  Raw per-read variant (before consensus)
    VC  Variant count in consensus call
    BC  Barcode sequence
    BQ  Mean barcode quality
    BK  Cluster centroid barcode
    CS  Cluster size
    VQ  Minimum variant quality (raw call)
        """)

    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("reference_fasta", help="Reference FASTA file")

    parser.add_argument("--promoter-region", required=True,
                        help="Promoter region as CHROM:START-END (1-based inclusive)")
    parser.add_argument("--barcode-region",
                        help="Barcode region START-END (1-based inclusive, same chrom as promoter)")
    parser.add_argument("--no-barcode", action="store_true", help="Skip barcode extraction")
    parser.add_argument("--snv-only", action="store_true", help="Exclude reads with indel variants in promoter")

    parser.add_argument("-o", "--output-dir", default="variant_calling", help="Output directory (default: variant_calling)")
    parser.add_argument("-n", "--sample-name", help="Sample name (default: from BAM filename)")
    parser.add_argument("--min-base-quality", type=int, default=20, help="Min base quality for variant calls (default: 20)")
    parser.add_argument("--save-pngs", action="store_true", help="Save individual PNG plots")
    parser.add_argument("--no-index", action="store_true", help="Skip BAM indexing")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress progress output")

    adv = parser.add_argument_group('Advanced options')
    adv.add_argument("--min-barcode-length", type=int, default=None, help="Min barcode length (default: region size - 1)")
    adv.add_argument("--max-barcode-length", type=int, default=None, help="Max barcode length (default: region size + 1)")
    adv.add_argument("--max-edit-distance", type=int, default=2, help="Max edit distance for barcode clustering (default: 2)")
    adv.add_argument("--min-jaccard", type=float, default=0.5, help="Min Jaccard index of variant sets for barcode merge (default: 0.5)")
    adv.add_argument("--variant-frequency-threshold", type=float, default=0.6, help="Min read fraction for variant consensus (default: 0.6)")
    adv.add_argument("--min-indel-reads", type=int, default=2, help="Min reads for indel consensus call (default: 2)")
    adv.add_argument("--top-variants-nucleosome", type=int, default=10, help="Top variants for nucleosome plots (default: 10)")
    adv.add_argument("--expected-variant-types", default="snv,ins,del", help="Variant types for coverage calc (default: snv,ins,del)")

    return parser.parse_args()


def main():
    args = parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.quiet:
        logging.getLogger().setLevel(logging.WARNING)

    promoter_chrom, promoter_start, promoter_end = parse_region_coords(args.promoter_region)
    if promoter_chrom is None:
        logger.error("--promoter-region must include chromosome: CHROM:START-END")
        sys.exit(1)
    extract_barcode_flag = not args.no_barcode
    barcode_start = barcode_end = None

    if extract_barcode_flag:
        if args.barcode_region is None:
            logger.error("--barcode-region is required unless --no-barcode is specified")
            sys.exit(1)
        _, barcode_start, barcode_end = parse_region_coords(args.barcode_region)

    if extract_barcode_flag:
        expected_bc_len = barcode_end - barcode_start + 1
        min_barcode_length = args.min_barcode_length if args.min_barcode_length is not None else max(expected_bc_len - 1, 1)
        max_barcode_length = args.max_barcode_length if args.max_barcode_length is not None else expected_bc_len + 1
    else:
        expected_bc_len = 15
        min_barcode_length = 14
        max_barcode_length = 16

    sample_name = args.sample_name if args.sample_name else Path(args.input_bam).stem
    for suffix in ['.sorted', '.aligned', '.merged', '_sorted', '_aligned', '_merged']:
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]

    output_dir = Path(args.output_dir)
    report_dir = output_dir / "report"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    output_bam = output_dir / f"{sample_name}.tagged.bam"

    # ===== PASS 1 =====
    try:
        stats, ref_seq, read_data, barcode_read_variants = process_bam_pass1(
            input_bam=args.input_bam, reference_fasta=args.reference_fasta,
            promoter_chrom=promoter_chrom, promoter_start=promoter_start, promoter_end=promoter_end,
            barcode_start=barcode_start, barcode_end=barcode_end, min_base_quality=args.min_base_quality,
            min_barcode_length=min_barcode_length, max_barcode_length=max_barcode_length,
            extract_barcode=extract_barcode_flag, snv_only=args.snv_only, quiet=args.quiet)
    except Exception as e:
        logger.error(f"Error processing BAM (pass 1): {e}")
        raise

    # ===== CLUSTERING + CONSENSUS =====
    cluster_resolutions = {}
    clustering_stats = None
    barcode_to_cluster = {}

    if extract_barcode_flag:
        clusters, barcode_to_cluster = cluster_barcodes(
            stats.barcode_variant_map, barcode_read_variants,
            max_edit_distance=args.max_edit_distance,
            expected_barcode_length=expected_bc_len,
            min_jaccard=args.min_jaccard,
            quiet=args.quiet)
        cluster_resolutions, clustering_stats = resolve_barcode_clusters(
            clusters, stats.barcode_variant_map, barcode_read_variants,
            variant_frequency_threshold=args.variant_frequency_threshold,
            min_indel_reads=args.min_indel_reads, quiet=args.quiet)

    # ===== PASS 2 =====
    cluster_corrections = {}
    if extract_barcode_flag:
        written, corrected, cluster_corrections = process_bam_pass2(
            input_bam=args.input_bam, output_bam=str(output_bam),
            promoter_chrom=promoter_chrom, promoter_start=promoter_start, promoter_end=promoter_end,
            barcode_start=barcode_start, barcode_end=barcode_end,
            read_data=read_data, barcode_to_cluster=barcode_to_cluster,
            cluster_resolutions=cluster_resolutions,
            extract_barcode=extract_barcode_flag, quiet=args.quiet)
    else:
        if not args.quiet:
            logger.info("Writing output BAM (no barcode)...")
        bam_in = pysam.AlignmentFile(args.input_bam, "rb")
        bam_out = pysam.AlignmentFile(str(output_bam), "wb", template=bam_in)
        promoter_start_0 = promoter_start - 1
        promoter_end_0 = promoter_end - 1
        written = 0
        try:
            for read in bam_in.fetch(promoter_chrom, promoter_start_0, promoter_end_0 + 1):
                if read.query_name in read_data:
                    rd = read_data[read.query_name]
                    read.set_tag("PV", rd['variant_tag'])
                    read.set_tag("VC", rd['num_variants'])
                    if rd['min_var_qual'] is not None:
                        read.set_tag("VQ", rd['min_var_qual'])
                    bam_out.write(read)
                    written += 1
        except ValueError as e:
            logger.warning(f"Could not fetch region: {e}")
        bam_in.close()
        bam_out.close()
        corrected = 0

    # ===== INDEX =====
    if not args.no_index:
        if not args.quiet:
            logger.info("Indexing output BAM...")
        try:
            pysam.index(str(output_bam))
        except Exception as e:
            logger.warning(f"Could not index BAM: {e}")

    # ===== REPORTS =====
    variant_types = [v.strip() for v in args.expected_variant_types.split(',')]
    promoter_length = promoter_end - promoter_start + 1
    coverage = calculate_expected_coverage(stats, promoter_length, ref_seq, variant_types)

    if not args.quiet:
        logger.info("Writing reports...")
    report_prefix = report_dir / sample_name
    write_stats_report(stats, f"{report_prefix}_stats.tsv", clustering_stats)
    write_variants_report(stats, f"{report_prefix}_variants.tsv")
    if extract_barcode_flag:
        write_cluster_resolution_report(cluster_resolutions, f"{report_prefix}_cluster_resolution.tsv")
    write_exclusion_report(stats, f"{report_prefix}_exclusion_report.tsv")
    write_expected_coverage_report(coverage, f"{report_prefix}_expected_coverage.tsv")
    if extract_barcode_flag and cluster_corrections:
        write_correction_report(cluster_corrections, f"{report_prefix}_corrections.tsv")

    if not args.quiet:
        logger.info("Generating plots...")
    figures = generate_plots(
        stats=stats, cluster_resolutions=cluster_resolutions, clustering_stats=clustering_stats,
        coverage=coverage, cluster_corrections=cluster_corrections,
        promoter_start=promoter_start, promoter_end=promoter_end,
        expected_barcode_length=expected_bc_len,
        top_variants_nucleosome=args.top_variants_nucleosome,
        output_prefix=str(report_prefix) if args.save_pngs else None,
        save_individual_pngs=args.save_pngs)

    generate_html_report(stats, cluster_resolutions, clustering_stats,
                          coverage, figures, f"{report_prefix}_report.html", sample_name,
                          include_barcode=extract_barcode_flag)
    generate_pdf_report(figures, f"{report_prefix}_report.pdf", sample_name)

    # ===== SUMMARY =====
    if not args.quiet:
        print("\n" + "=" * 60)
        print("MPRA VARIANT-BARCODE TAGGER SUMMARY")
        print("=" * 60)
        print(f"Total reads processed:     {stats.total_reads:>12,}")
        if stats.total_reads > 0:
            print(f"Reads passing filters:     {stats.passed_reads:>12,} ({stats.passed_reads/stats.total_reads*100:.1f}%)")
        if args.snv_only:
            print(f"  (SNV-only: {stats.exclusion_counts.get(ExclusionReason.HAS_INDEL_VARIANT, 0):,} reads excluded for indels)")
        print(f"  - Wild-type reads:       {stats.wt_reads:>12,}")
        print(f"  - Single variant reads:  {stats.single_variant_reads:>12,}")
        print(f"  - Multi-variant reads:   {stats.multi_variant_reads:>12,}")
        print(f"Unique variants found:     {len(stats.variant_counts):>12,}")

        if extract_barcode_flag and clustering_stats:
            print(f"\nUnique barcodes found:     {len(stats.barcode_variant_map):>12,}")
            print(f"\n--- Barcode Clustering ---")
            print(f"Clusters after merging:    {clustering_stats.total_clusters_after:>12,}")
            print(f"Barcodes merged:           {clustering_stats.barcodes_merged:>12,}")
            print(f"Singleton clusters:        {clustering_stats.singleton_clusters:>12,}")
            print(f"  - Resolved clusters:     {clustering_stats.resolved_clusters:>12,}")
            print(f"  - Ambiguous clusters:    {clustering_stats.ambiguous_clusters:>12,}")
            print(f"  - Low-confidence:        {clustering_stats.low_confidence_clusters:>12,}")
            print(f"Reads with corrected tags: {corrected:>12,}")
            if cluster_corrections:
                print(f"  Clusters with corrections: {len(cluster_corrections):>8,}")
                all_removed = sum(len(cc.variants_removed) for cc in cluster_corrections.values())
                print(f"  Total variants removed:    {all_removed:>8,}")

        print("-" * 60)
        print("Expected Variant Coverage:")
        print(f"  - SNVs:                  {coverage.observed_snvs:>6,} / {coverage.expected_snvs:>6,} ({coverage.snv_coverage:.1f}%)")
        print(f"  - Insertions (1bp):      {coverage.observed_insertions:>6,} / {coverage.expected_insertions:>6,} ({coverage.insertion_coverage:.1f}%)")
        print(f"  - Deletions (1bp):       {coverage.observed_deletions:>6,} / {coverage.expected_deletions:>6,} ({coverage.deletion_coverage:.1f}%)")
        print("=" * 60)
        print(f"\nOutput directory: {output_dir}")
        print(f"  - Tagged BAM: {output_bam.name}")
        print(f"  - Reports: {report_dir.name}/")

    logger.info("Done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
