#!/usr/bin/env python3
"""
MPRA Variant-Barcode Tagger

Tags PacBio MPRA reads with promoter variants and barcode sequences.
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


class ExclusionReason(Enum):
    PASS = "pass"
    UNMAPPED = "unmapped"
    SECONDARY = "secondary_alignment"
    SUPPLEMENTARY = "supplementary_alignment"
    NO_PROMOTER_COVERAGE = "no_promoter_coverage"
    NO_BARCODE_COVERAGE = "no_barcode_coverage"
    BARCODE_TOO_SHORT = "barcode_too_short"
    LOW_VARIANT_QUALITY = "low_variant_quality"
    FORBIDDEN_CIGAR = "forbidden_cigar_operation"


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
    variant_qualities: list = field(default_factory=list)
    barcode_qualities: list = field(default_factory=list)
    variants_per_read: list = field(default_factory=list)
    variant_nucleosome_counts: dict = field(default_factory=lambda: defaultdict(list))
    snv_position_counts: Counter = field(default_factory=Counter)
    ins_position_counts: Counter = field(default_factory=Counter)
    del_position_counts: Counter = field(default_factory=Counter)
    observed_snvs: set = field(default_factory=set)
    observed_insertions: set = field(default_factory=set)
    observed_deletions: set = field(default_factory=set)
    single_variant_snvs: set = field(default_factory=set)  # SNVs from reads with only one variant


@dataclass
class BarcodeResolution:
    consensus_variant: str
    is_ambiguous: bool
    vote_proportion: float
    total_reads: int
    competing_variants: dict


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


def parse_promoter_variants(read, ref_seq, region_start_0based, region_end_0based, min_quality=0):
    """Parse CIGAR string to extract variants in the promoter region."""
    variants = []
    
    if read.reference_end is None:
        return [], ExclusionReason.NO_PROMOTER_COVERAGE
    if read.reference_start > region_end_0based or read.reference_end <= region_start_0based:
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
        
        if op == 0:  # M - forbidden in PacBio
            return [], ExclusionReason.FORBIDDEN_CIGAR
        
        elif op == 7:  # = (sequence match)
            query_idx += length
            ref_pos += length
        
        elif op == 8:  # X (sequence mismatch) - SNVs
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
                            variants.append(VariantCall(
                                position=abs_ref_pos + 1, ref=ref_base, alt=alt_base,
                                var_type=VariantType.SNV, quality=qual
                            ))
            query_idx += length
            ref_pos += length
        
        elif op == 1:  # I (insertion)
            if region_start_0based <= ref_pos <= region_end_0based:
                inserted_seq = query_seq[query_idx:query_idx + length]
                if query_quals is not None:
                    ins_quals = query_quals[query_idx:query_idx + length]
                    min_ins_qual = min(ins_quals) if ins_quals else None
                else:
                    min_ins_qual = None
                
                if min_ins_qual is None or min_ins_qual >= min_quality:
                    variants.append(VariantCall(
                        position=ref_pos + 1, ref="", alt=inserted_seq,
                        var_type=VariantType.INSERTION, quality=min_ins_qual
                    ))
            query_idx += length
        
        elif op == 2:  # D (deletion)
            if ref_pos + length > region_start_0based and ref_pos <= region_end_0based:
                overlap_start = max(ref_pos, region_start_0based)
                overlap_end = min(ref_pos + length, region_end_0based + 1)
                
                if overlap_start < overlap_end:
                    rel_start = overlap_start - region_start_0based
                    rel_end = overlap_end - region_start_0based
                    deleted_seq = ref_seq[rel_start:rel_end]
                    variants.append(VariantCall(
                        position=overlap_start + 1, ref=deleted_seq, alt="",
                        var_type=VariantType.DELETION, quality=None
                    ))
            ref_pos += length
        
        elif op == 4:  # S (soft clip)
            query_idx += length
        elif op == 5:  # H (hard clip)
            pass
        elif op == 3:  # N (reference skip)
            ref_pos += length
    
    return variants, ExclusionReason.PASS


def extract_barcode_seq(read, ref_start_0based, ref_end_0based, expected_length=15, min_length=13):
    """Extract barcode sequence from read at specified reference coordinates."""
    if read.is_unmapped:
        return None, None, ExclusionReason.UNMAPPED
    if read.reference_end is None:
        return None, None, ExclusionReason.NO_BARCODE_COVERAGE
    if read.reference_start > ref_end_0based or read.reference_end <= ref_start_0based:
        return None, None, ExclusionReason.NO_BARCODE_COVERAGE
    
    aligned_pairs = read.get_aligned_pairs(matches_only=False)
    barcode_query_positions = []
    in_region = False
    region_started = False
    
    for query_pos, ref_pos in aligned_pairs:
        if ref_pos is not None and ref_start_0based <= ref_pos < ref_end_0based:
            in_region = True
            region_started = True
            if query_pos is not None:
                barcode_query_positions.append(query_pos)
        elif in_region and ref_pos is None and query_pos is not None:
            barcode_query_positions.append(query_pos)
        elif region_started and ref_pos is not None and ref_pos >= ref_end_0based:
            break
        elif not region_started and ref_pos is not None and ref_pos >= ref_end_0based:
            return None, None, ExclusionReason.NO_BARCODE_COVERAGE
    
    if len(barcode_query_positions) < min_length:
        if len(barcode_query_positions) == 0:
            return None, None, ExclusionReason.NO_BARCODE_COVERAGE
        else:
            return None, None, ExclusionReason.BARCODE_TOO_SHORT
    
    barcode_seq = ''.join([read.query_sequence[pos] for pos in barcode_query_positions])
    
    if read.query_qualities is not None:
        barcode_quals = [read.query_qualities[pos] for pos in barcode_query_positions]
        mean_quality = sum(barcode_quals) / len(barcode_quals)
    else:
        mean_quality = None
    
    return barcode_seq, mean_quality, ExclusionReason.PASS


def resolve_barcode_variants(barcode_variant_map, consensus_threshold=0.9):
    """Resolve barcode-to-variant mapping using majority voting."""
    resolutions = {}
    for barcode, variant_counts in barcode_variant_map.items():
        total_reads = sum(variant_counts.values())
        if total_reads == 0:
            continue
        most_common = variant_counts.most_common(1)[0]
        consensus_variant, consensus_count = most_common
        vote_proportion = consensus_count / total_reads
        is_ambiguous = vote_proportion < consensus_threshold
        resolutions[barcode] = BarcodeResolution(
            consensus_variant=consensus_variant, is_ambiguous=is_ambiguous,
            vote_proportion=vote_proportion, total_reads=total_reads,
            competing_variants=dict(variant_counts)
        )
    return resolutions


def variants_to_tag(variants):
    """Convert list of variants to tag string."""
    if not variants:
        return "WT"
    variant_ids = sorted([v.id for v in variants])
    return json.dumps(variant_ids)


def calculate_expected_coverage(stats, promoter_length, ref_seq, variant_types=None):
    """Calculate expected vs observed variant coverage."""
    if variant_types is None:
        variant_types = ['snv', 'ins', 'del']
    
    expected_snvs = promoter_length * 3 if 'snv' in variant_types else 0
    expected_insertions = promoter_length * 4 if 'ins' in variant_types else 0
    expected_deletions = promoter_length if 'del' in variant_types else 0
    
    observed_snvs = len(stats.observed_snvs) if 'snv' in variant_types else 0
    
    observed_insertions = 0
    if 'ins' in variant_types:
        for var_id in stats.observed_insertions:
            if '+' in var_id:
                seq = var_id.split('+')[1]
                if len(seq) == 1:
                    observed_insertions += 1
    
    observed_deletions = 0
    if 'del' in variant_types:
        for var_id in stats.observed_deletions:
            parts = var_id.split(':')
            if len(parts) == 2:
                length_str = ''.join(c for c in parts[1] if c.isdigit())
                if length_str and int(length_str) == 1:
                    observed_deletions += 1
    
    snv_coverage = (observed_snvs / expected_snvs * 100) if expected_snvs > 0 else 0
    insertion_coverage = (observed_insertions / expected_insertions * 100) if expected_insertions > 0 else 0
    deletion_coverage = (observed_deletions / expected_deletions * 100) if expected_deletions > 0 else 0
    
    return ExpectedCoverage(
        promoter_length=promoter_length,
        expected_snvs=expected_snvs, observed_snvs=observed_snvs, snv_coverage=snv_coverage,
        expected_insertions=expected_insertions, observed_insertions=observed_insertions, insertion_coverage=insertion_coverage,
        expected_deletions=expected_deletions, observed_deletions=observed_deletions, deletion_coverage=deletion_coverage
    )


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
        if op == 0:  # M - forbidden in PacBio
            return [], ExclusionReason.FORBIDDEN_CIGAR
        elif op == 7:  # = (match)
            query_idx += length
            ref_pos += length
        elif op == 8:  # X (mismatch/SNV)
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
        elif op == 1:  # I (insertion)
            if region_start_0based <= ref_pos <= region_end_0based:
                inserted_seq = query_seq[query_idx:query_idx + length]
                min_ins_qual = min(query_quals[query_idx:query_idx + length]) if query_quals is not None else None
                if min_ins_qual is None or min_ins_qual >= min_quality:
                    variants.append(VariantCall(position=ref_pos + 1, ref="", alt=inserted_seq, var_type=VariantType.INSERTION, quality=min_ins_qual))
            query_idx += length
        elif op == 2:  # D (deletion)
            if ref_pos + length > region_start_0based and ref_pos <= region_end_0based:
                overlap_start = max(ref_pos, region_start_0based)
                overlap_end = min(ref_pos + length, region_end_0based + 1)
                if overlap_start < overlap_end:
                    rel_start = overlap_start - region_start_0based
                    rel_end = overlap_end - region_start_0based
                    deleted_seq = ref_seq[rel_start:rel_end]
                    variants.append(VariantCall(position=overlap_start + 1, ref=deleted_seq, alt="", var_type=VariantType.DELETION, quality=None))
            ref_pos += length
        elif op == 4:  # S (soft clip)
            query_idx += length
        elif op == 5:  # H (hard clip)
            pass
        elif op == 3:  # N (ref skip)
            ref_pos += length
    return variants, ExclusionReason.PASS


def extract_barcode_seq(read, ref_start_0based, ref_end_0based, expected_length=15, min_length=13):
    """Extract barcode sequence from read at specified reference coordinates."""
    if read.is_unmapped:
        return None, None, ExclusionReason.UNMAPPED
    if read.reference_end is None or read.reference_start > ref_end_0based or read.reference_end <= ref_start_0based:
        return None, None, ExclusionReason.NO_BARCODE_COVERAGE
    aligned_pairs = read.get_aligned_pairs(matches_only=False)
    barcode_query_positions = []
    in_region = False
    region_started = False
    for query_pos, ref_pos in aligned_pairs:
        if ref_pos is not None and ref_start_0based <= ref_pos < ref_end_0based:
            in_region = True
            region_started = True
            if query_pos is not None:
                barcode_query_positions.append(query_pos)
        elif in_region and ref_pos is None and query_pos is not None:
            barcode_query_positions.append(query_pos)
        elif region_started and ref_pos is not None and ref_pos >= ref_end_0based:
            break
        elif not region_started and ref_pos is not None and ref_pos >= ref_end_0based:
            return None, None, ExclusionReason.NO_BARCODE_COVERAGE
    if len(barcode_query_positions) < min_length:
        return (None, None, ExclusionReason.NO_BARCODE_COVERAGE) if len(barcode_query_positions) == 0 else (None, None, ExclusionReason.BARCODE_TOO_SHORT)
    barcode_seq = ''.join([read.query_sequence[pos] for pos in barcode_query_positions])
    mean_quality = sum(read.query_qualities[pos] for pos in barcode_query_positions) / len(barcode_query_positions) if read.query_qualities is not None else None
    return barcode_seq, mean_quality, ExclusionReason.PASS


def resolve_barcode_variants(barcode_variant_map, consensus_threshold=0.9):
    """Resolve barcode-to-variant mapping using majority voting."""
    resolutions = {}
    for barcode, variant_counts in barcode_variant_map.items():
        total_reads = sum(variant_counts.values())
        if total_reads == 0:
            continue
        most_common = variant_counts.most_common(1)[0]
        consensus_variant, consensus_count = most_common
        vote_proportion = consensus_count / total_reads
        resolutions[barcode] = BarcodeResolution(consensus_variant=consensus_variant, is_ambiguous=vote_proportion < consensus_threshold, vote_proportion=vote_proportion, total_reads=total_reads, competing_variants=dict(variant_counts))
    return resolutions


def variants_to_tag(variants):
    if not variants:
        return "WT"
    return json.dumps(sorted([v.id for v in variants]))


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


def process_bam(input_bam, output_bam, reference_fasta, promoter_chrom, promoter_start, promoter_end,
                barcode_start=None, barcode_end=None, min_base_quality=20, expected_barcode_length=15, 
                min_barcode_length=13, extract_barcode=True, quiet=False):
    """Main BAM processing function."""
    stats = ProcessingStats()
    promoter_start_0 = promoter_start - 1
    promoter_end_0 = promoter_end - 1
    
    # Handle barcode coordinates
    if extract_barcode:
        if barcode_start is None or barcode_end is None:
            raise ValueError("barcode_start and barcode_end are required when extract_barcode=True")
        barcode_start_0 = barcode_start - 1
        barcode_end_0 = barcode_end
    else:
        barcode_start_0 = None
        barcode_end_0 = None
    
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)
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
        logger.info(f"Processing BAM: {input_bam}")
        logger.info(f"Promoter region: {promoter_chrom}:{promoter_start}-{promoter_end} ({promoter_end - promoter_start + 1} bp)")
        if extract_barcode:
            logger.info(f"Barcode region: {promoter_chrom}:{barcode_start}-{barcode_end}")
        else:
            logger.info("Barcode extraction: disabled")
    
    # Determine fetch region
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
            
            variants, var_status = parse_promoter_variants(read, ref_seq, promoter_start_0, promoter_end_0, min_base_quality)
            if var_status != ExclusionReason.PASS:
                stats.exclusion_counts[var_status] += 1
                continue
            
            # Extract barcode if enabled
            barcode = None
            barcode_qual = None
            if extract_barcode:
                barcode, barcode_qual, bc_status = extract_barcode_seq(read, barcode_start_0, barcode_end_0, expected_barcode_length, min_barcode_length)
                if bc_status != ExclusionReason.PASS:
                    stats.exclusion_counts[bc_status] += 1
                    continue
            
            nucleosome_count = read.get_tag("nc") if read.has_tag("nc") else None
            variant_tag = variants_to_tag(variants)
            stats.passed_reads += 1
            num_variants = len(variants)
            stats.variants_per_read.append(num_variants)
            
            if num_variants == 0:
                stats.wt_reads += 1
            elif num_variants == 1:
                stats.single_variant_reads += 1
                # Track SNVs that appear as the only variant in a read
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
                if var.quality is not None:
                    stats.variant_qualities.append(var.quality)
            
            # Track barcode stats only if barcode extraction is enabled
            if extract_barcode and barcode:
                stats.barcode_lengths[len(barcode)] += 1
                if barcode_qual is not None:
                    stats.barcode_qualities.append(barcode_qual)
                stats.barcode_variant_map[barcode][variant_tag] += 1
            
            if nucleosome_count is not None:
                stats.variant_nucleosome_counts[variant_tag].append(nucleosome_count)
            
            read.set_tag("PV", variant_tag)
            read.set_tag("VC", num_variants)
            if extract_barcode and barcode:
                read.set_tag("BC", barcode)
                if barcode_qual is not None:
                    read.set_tag("BQ", int(round(barcode_qual)))
            if variants:
                min_var_qual = min((v.quality for v in variants if v.quality is not None), default=None)
                if min_var_qual is not None:
                    read.set_tag("VQ", min_var_qual)
            bam_out.write(read)
            
            if not quiet and stats.total_reads % 100000 == 0:
                logger.info(f"Processed {stats.total_reads:,} reads, {stats.passed_reads:,} passed...")
    except ValueError as e:
        logger.warning(f"Could not fetch region: {e}")
    
    bam_in.close()
    bam_out.close()
    fasta.close()
    if not quiet:
        logger.info(f"Finished processing {stats.total_reads:,} reads, wrote {stats.passed_reads:,} to output BAM")
    return stats, ref_seq


def write_stats_report(stats, output_path):
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
        f.write(f"total_snvs\t{stats.variant_type_counts.get('snv', 0)}\n")
        f.write(f"total_insertions\t{stats.variant_type_counts.get('insertion', 0)}\n")
        f.write(f"total_deletions\t{stats.variant_type_counts.get('deletion', 0)}\n")
        f.write(f"transitions\t{stats.transition_count}\n")
        f.write(f"transversions\t{stats.transversion_count}\n")
        ti_tv = stats.transition_count / stats.transversion_count if stats.transversion_count > 0 else float('inf')
        f.write(f"ti_tv_ratio\t{ti_tv:.4f}\n")
        f.write(f"unique_variants\t{len(stats.variant_counts)}\n")
        f.write(f"unique_barcodes\t{len(stats.barcode_variant_map)}\n")


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


def write_barcode_variant_map(stats, barcode_resolution, output_path):
    with open(output_path, 'w') as f:
        f.write("barcode\tconsensus_variant\tis_ambiguous\tvote_proportion\ttotal_reads\tall_variants\n")
        for barcode in sorted(stats.barcode_variant_map.keys()):
            if barcode in barcode_resolution:
                res = barcode_resolution[barcode]
                f.write(f"{barcode}\t{res.consensus_variant}\t{res.is_ambiguous}\t{res.vote_proportion:.4f}\t{res.total_reads}\t{json.dumps(res.competing_variants)}\n")


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


def generate_plots(stats, barcode_resolution, coverage, promoter_start, promoter_end, 
                   top_variants_nucleosome=10, output_prefix=None, save_individual_pngs=False):
    """Generate diagnostic plots. Returns dict of figure objects."""
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
    
    # 1. Variants per read
    if stats.variants_per_read:
        fig, ax = plt.subplots()
        max_var = max(stats.variants_per_read)
        ax.hist(stats.variants_per_read, bins=range(0, min(max_var + 2, 20)), edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xlabel('Number of Variants per Read')
        ax.set_ylabel('Read Count')
        ax.set_title('Distribution of Variants per Read')
        total = len(stats.variants_per_read)
        ax.text(0.95, 0.95, f'WT: {stats.wt_reads/total*100:.1f}%\nSingle: {stats.single_variant_reads/total*100:.1f}%',
                transform=ax.transAxes, ha='right', va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        figures['variants_per_read'] = fig
        if save_individual_pngs and output_prefix:
            plt.savefig(f"{output_prefix}_variants_per_read.png", dpi=150)
    
    # 2. Barcodes per variant
    if stats.barcode_variant_map:
        variant_to_barcodes = defaultdict(set)
        for barcode, var_counts in stats.barcode_variant_map.items():
            for var_tag in var_counts.keys():
                variant_to_barcodes[var_tag].add(barcode)
        bc_counts = [len(bcs) for bcs in variant_to_barcodes.values()]
        if bc_counts:
            fig, ax = plt.subplots()
            ax.hist(bc_counts, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
            ax.set_xlabel('Number of Unique Barcodes')
            ax.set_ylabel('Variant Count')
            ax.set_title('Distribution of Barcodes per Variant')
            ax.set_yscale('log')
            ax.axvline(np.mean(bc_counts), color='red', linestyle='--', label=f'Mean: {np.mean(bc_counts):.1f}')
            ax.legend()
            plt.tight_layout()
            figures['barcodes_per_variant'] = fig
            if save_individual_pngs and output_prefix:
                plt.savefig(f"{output_prefix}_barcodes_per_variant.png", dpi=150)
    
    # 3. Barcode resolution
    if barcode_resolution:
        resolved = sum(1 for r in barcode_resolution.values() if not r.is_ambiguous)
        ambiguous = sum(1 for r in barcode_resolution.values() if r.is_ambiguous)
        fig, ax = plt.subplots()
        ax.bar(['Resolved', 'Ambiguous'], [resolved, ambiguous], color=['forestgreen', 'firebrick'], alpha=0.7, edgecolor='black')
        ax.set_ylabel('Number of Barcodes')
        ax.set_title('Barcode-Variant Resolution Status')
        plt.tight_layout()
        figures['barcode_resolution'] = fig
        if save_individual_pngs and output_prefix:
            plt.savefig(f"{output_prefix}_barcode_resolution.png", dpi=150)
    
    # 4. SNV read counts
    snv_variants = {k: v for k, v in stats.variant_counts.items() if '>' in k}
    if snv_variants:
        fig, ax = plt.subplots(figsize=(12, 6))
        sorted_snvs = sorted(snv_variants.items(), key=lambda x: -x[1])[:50]
        ax.bar(range(len(sorted_snvs)), [v[1] for v in sorted_snvs], edgecolor='black', alpha=0.7, color='steelblue')
        ax.set_xlabel('SNV Variant')
        ax.set_ylabel('Read Count')
        ax.set_title(f'Reads per SNV Variant (Top {len(sorted_snvs)})')
        ax.set_xticks(range(len(sorted_snvs)))
        ax.set_xticklabels([v[0] for v in sorted_snvs], rotation=90, fontsize=8)
        plt.tight_layout()
        figures['snv_read_counts'] = fig
        if save_individual_pngs and output_prefix:
            plt.savefig(f"{output_prefix}_snv_read_counts.png", dpi=150)
    
    # 5. Insertion read counts
    ins_variants = {k: v for k, v in stats.variant_counts.items() if '+' in k}
    if ins_variants:
        fig, ax = plt.subplots(figsize=(12, 6))
        sorted_ins = sorted(ins_variants.items(), key=lambda x: -x[1])[:50]
        ax.bar(range(len(sorted_ins)), [v[1] for v in sorted_ins], edgecolor='black', alpha=0.7, color='darkorange')
        ax.set_xlabel('Insertion Variant')
        ax.set_ylabel('Read Count')
        ax.set_title(f'Reads per Insertion Variant (Top {len(sorted_ins)})')
        ax.set_xticks(range(len(sorted_ins)))
        ax.set_xticklabels([v[0] for v in sorted_ins], rotation=90, fontsize=8)
        plt.tight_layout()
        figures['insertion_read_counts'] = fig
        if save_individual_pngs and output_prefix:
            plt.savefig(f"{output_prefix}_insertion_read_counts.png", dpi=150)
    
    # 6. Deletion read counts
    del_variants = {k: v for k, v in stats.variant_counts.items() if '>' not in k and '+' not in k}
    if del_variants:
        fig, ax = plt.subplots(figsize=(12, 6))
        sorted_del = sorted(del_variants.items(), key=lambda x: -x[1])[:50]
        ax.bar(range(len(sorted_del)), [v[1] for v in sorted_del], edgecolor='black', alpha=0.7, color='firebrick')
        ax.set_xlabel('Deletion Variant')
        ax.set_ylabel('Read Count')
        ax.set_title(f'Reads per Deletion Variant (Top {len(sorted_del)})')
        ax.set_xticks(range(len(sorted_del)))
        ax.set_xticklabels([v[0] for v in sorted_del], rotation=90, fontsize=8)
        plt.tight_layout()
        figures['deletion_read_counts'] = fig
        if save_individual_pngs and output_prefix:
            plt.savefig(f"{output_prefix}_deletion_read_counts.png", dpi=150)
    
    # 7. Position coverage
    if stats.snv_position_counts or stats.ins_position_counts or stats.del_position_counts:
        fig, ax = plt.subplots(figsize=(14, 6))
        positions = list(range(promoter_start, promoter_end + 1))
        snv_counts = [stats.snv_position_counts.get(p, 0) for p in positions]
        ins_counts = [stats.ins_position_counts.get(p, 0) for p in positions]
        del_counts = [stats.del_position_counts.get(p, 0) for p in positions]
        ax.bar(positions, snv_counts, label='SNVs', alpha=0.7, color='steelblue')
        ax.bar(positions, ins_counts, bottom=snv_counts, label='Insertions', alpha=0.7, color='darkorange')
        ax.bar(positions, del_counts, bottom=[s+i for s,i in zip(snv_counts, ins_counts)], label='Deletions', alpha=0.7, color='firebrick')
        ax.set_xlabel('Position in Promoter')
        ax.set_ylabel('Variant Count')
        ax.set_title('Variant Frequency by Position')
        ax.legend()
        plt.tight_layout()
        figures['position_coverage'] = fig
        if save_individual_pngs and output_prefix:
            plt.savefig(f"{output_prefix}_position_coverage.png", dpi=150)
    
    # 8. Position frequency
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
        if save_individual_pngs and output_prefix:
            plt.savefig(f"{output_prefix}_position_frequency.png", dpi=150)
    
    # 9. Expected coverage
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
    if save_individual_pngs and output_prefix:
        plt.savefig(f"{output_prefix}_expected_coverage.png", dpi=150)
    
    # 10. SNV coverage comparison: all reads vs single-variant reads only
    if coverage.expected_snvs > 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        categories = ['All Reads', 'Single-Variant\nReads Only']
        all_snvs = len(stats.observed_snvs)
        single_var_snvs = len(stats.single_variant_snvs)
        expected_snvs = coverage.expected_snvs
        
        x = np.arange(len(categories))
        width = 0.5
        bars = ax.bar(x, [all_snvs, single_var_snvs], width, alpha=0.7, color=['steelblue', 'forestgreen'], edgecolor='black')
        ax.axhline(expected_snvs, color='red', linestyle='--', linewidth=2, label=f'Expected: {expected_snvs:,}')
        ax.set_ylabel('Unique SNVs Observed')
        ax.set_title('SNV Coverage: All Reads vs Single-Variant Reads')
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.legend()
        
        # Add percentage labels
        for i, (count, bar) in enumerate(zip([all_snvs, single_var_snvs], bars)):
            pct = count / expected_snvs * 100
            ax.annotate(f'{count:,}\n({pct:.1f}%)', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                       ha='center', va='bottom', fontsize=11, fontweight='bold')
        
        plt.tight_layout()
        figures['snv_coverage_comparison'] = fig
        if save_individual_pngs and output_prefix:
            plt.savefig(f"{output_prefix}_snv_coverage_comparison.png", dpi=150)
    
    # 10. Nucleosome distribution
    if stats.variant_nucleosome_counts:
        variant_read_counts = {tag: len(counts) for tag, counts in stats.variant_nucleosome_counts.items()}
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
                    nc_counts = stats.variant_nucleosome_counts[var_tag]
                    if nc_counts:
                        ax.hist(nc_counts, bins=range(0, max(nc_counts) + 2), edgecolor='black', alpha=0.7, color='steelblue')
                        ax.set_xlabel('Nucleosome Count')
                        ax.set_ylabel('Read Count')
                        label = var_tag if len(var_tag) < 30 else var_tag[:27] + '...'
                        ax.set_title(f'{label}\n(n={len(nc_counts)})')
                        ax.axvline(np.mean(nc_counts), color='red', linestyle='--', label=f'Mean: {np.mean(nc_counts):.1f}')
                        ax.legend(fontsize=8)
            for idx in range(len(variants_to_plot), len(axes)):
                axes[idx].set_visible(False)
            plt.suptitle('Nucleosome Count Distribution by Variant', fontsize=14, y=1.02)
            plt.tight_layout()
            figures['nucleosome_distribution'] = fig
            if save_individual_pngs and output_prefix:
                plt.savefig(f"{output_prefix}_nucleosome_distribution.png", dpi=150, bbox_inches='tight')
    
    return figures


def fig_to_base64(fig):
    """Convert matplotlib figure to base64 string."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def generate_html_report(stats, barcode_resolution, coverage, figures, output_path, sample_name, include_barcode=True):
    """Generate HTML report with embedded plots."""
    pass_rate = stats.passed_reads / stats.total_reads * 100 if stats.total_reads > 0 else 0
    resolved = sum(1 for r in barcode_resolution.values() if not r.is_ambiguous) if include_barcode else 0
    ambiguous = sum(1 for r in barcode_resolution.values() if r.is_ambiguous) if include_barcode else 0
    ti_tv = stats.transition_count / stats.transversion_count if stats.transversion_count > 0 else float('inf')
    wt_pct = stats.wt_reads/stats.passed_reads*100 if stats.passed_reads > 0 else 0
    
    # Build summary cards based on whether barcode is included
    barcode_card = ""
    if include_barcode:
        barcode_card = f'<div class="card"><h3>Unique Barcodes</h3><div class="value">{len(stats.barcode_variant_map):,}</div><div class="subvalue">{resolved:,} resolved, {ambiguous:,} ambiguous</div></div>'
    
    report_title = "MPRA Variant-Barcode QC Report" if include_barcode else "MPRA Variant QC Report"
    
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
.plot-container img {{ max-width: 100%; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
.plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }}
table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
th {{ background-color: #f8f9fa; font-weight: 600; }}
.good {{ color: #28a745; }} .warning {{ color: #ffc107; }} .bad {{ color: #dc3545; }}
</style></head><body>
<div class="header"><h1>{report_title}</h1><p>Sample: {sample_name}</p></div>
<div class="summary-cards">
<div class="card"><h3>Total Reads</h3><div class="value">{stats.total_reads:,}</div></div>
<div class="card"><h3>Passed Reads</h3><div class="value">{stats.passed_reads:,}</div><div class="subvalue">{pass_rate:.1f}% pass rate</div></div>
<div class="card"><h3>Unique Variants</h3><div class="value">{len(stats.variant_counts):,}</div></div>
{barcode_card}
<div class="card"><h3>WT Reads</h3><div class="value">{stats.wt_reads:,}</div><div class="subvalue">{wt_pct:.1f}% of passed</div></div>
<div class="card"><h3>Ti/Tv Ratio</h3><div class="value">{ti_tv:.2f}</div><div class="subvalue">{stats.transition_count:,} Ti / {stats.transversion_count:,} Tv</div></div>
</div>
<div class="section"><h2>Read Classification</h2>
<table><tr><th>Category</th><th>Count</th><th>Percentage</th></tr>
<tr><td>Wild-type reads</td><td>{stats.wt_reads:,}</td><td>{stats.wt_reads/stats.passed_reads*100:.2f}%</td></tr>
<tr><td>Single variant reads</td><td>{stats.single_variant_reads:,}</td><td>{stats.single_variant_reads/stats.passed_reads*100:.2f}%</td></tr>
<tr><td>Multi-variant reads</td><td>{stats.multi_variant_reads:,}</td><td>{stats.multi_variant_reads/stats.passed_reads*100:.2f}%</td></tr>
</table><div class="plot-grid">"""
    
    for plot_name in ['variants_per_read', 'read_classification']:
        if plot_name in figures:
            html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures[plot_name])}" alt="{plot_name}"></div>'
    
    html += f"""</div></div>
<div class="section"><h2>Expected Variant Coverage</h2>
<table><tr><th>Variant Type</th><th>Expected</th><th>Observed</th><th>Coverage</th></tr>
<tr><td>SNVs</td><td>{coverage.expected_snvs:,}</td><td>{coverage.observed_snvs:,}</td>
<td class="{'good' if coverage.snv_coverage > 80 else 'warning' if coverage.snv_coverage > 50 else 'bad'}">{coverage.snv_coverage:.1f}%</td></tr>
<tr><td>Single-bp Insertions</td><td>{coverage.expected_insertions:,}</td><td>{coverage.observed_insertions:,}</td>
<td class="{'good' if coverage.insertion_coverage > 80 else 'warning' if coverage.insertion_coverage > 50 else 'bad'}">{coverage.insertion_coverage:.1f}%</td></tr>
<tr><td>Single-bp Deletions</td><td>{coverage.expected_deletions:,}</td><td>{coverage.observed_deletions:,}</td>
<td class="{'good' if coverage.deletion_coverage > 80 else 'warning' if coverage.deletion_coverage > 50 else 'bad'}">{coverage.deletion_coverage:.1f}%</td></tr>
</table>"""
    
    if 'expected_coverage' in figures:
        html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures["expected_coverage"])}" alt="expected_coverage"></div>'
    
    if 'snv_coverage_comparison' in figures:
        html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures["snv_coverage_comparison"])}" alt="snv_coverage_comparison"></div>'
    
    html += '</div><div class="section"><h2>Variant Distribution</h2><div class="plot-grid">'
    for plot_name in ['position_coverage', 'position_frequency']:
        if plot_name in figures:
            html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures[plot_name])}" alt="{plot_name}"></div>'
    
    html += '</div></div><div class="section"><h2>Reads per Variant Type</h2><div class="plot-grid">'
    for plot_name in ['snv_read_counts', 'insertion_read_counts', 'deletion_read_counts']:
        if plot_name in figures:
            html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures[plot_name])}" alt="{plot_name}"></div>'
    
    html += '</div></div>'
    
    # Only include barcode analysis section if barcode extraction was enabled
    if include_barcode and ('barcodes_per_variant' in figures or 'barcode_resolution' in figures):
        html += '<div class="section"><h2>Barcode Analysis</h2><div class="plot-grid">'
        for plot_name in ['barcodes_per_variant', 'barcode_resolution']:
            if plot_name in figures:
                html += f'<div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures[plot_name])}" alt="{plot_name}"></div>'
        html += '</div></div>'
    
    if 'nucleosome_distribution' in figures:
        html += f'<div class="section"><h2>Nucleosome Distribution by Variant</h2><div class="plot-container"><img src="data:image/png;base64,{fig_to_base64(figures["nucleosome_distribution"])}" alt="nucleosome_distribution"></div></div>'
    
    html += '</body></html>'
    
    with open(output_path, 'w') as f:
        f.write(html)
    logger.info(f"Generated HTML report: {output_path}")


def generate_pdf_report(figures, output_path, sample_name):
    """Generate PDF report with all plots."""
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


def parse_args():
    parser = argparse.ArgumentParser(description="Tag MPRA reads with promoter variants and barcodes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python mpra_variant_barcode_tagger.py input.bam reference.fasta \\
        --promoter-chrom plasmid --promoter-start 100 --promoter-end 500 \\
        --barcode-start 3294 --barcode-end 3309

Output structure:
    variant_calling/
    ├── {sample}.tagged.bam
    ├── {sample}.tagged.bam.bai
    └── report/
        ├── {sample}_stats.tsv
        ├── {sample}_variants.tsv
        ├── {sample}_report.html
        └── {sample}_report.pdf
        """)
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("reference_fasta", help="Reference FASTA file")
    parser.add_argument("--promoter-chrom", required=True, help="Chromosome/contig name for promoter region")
    parser.add_argument("--promoter-start", type=int, required=True, help="Start position of promoter region (1-based)")
    parser.add_argument("--promoter-end", type=int, required=True, help="End position of promoter region (1-based)")
    parser.add_argument("--barcode-start", type=int, help="Start position of barcode region (1-based). Required unless --no-barcode is set.")
    parser.add_argument("--barcode-end", type=int, help="End position of barcode region (1-based). Required unless --no-barcode is set.")
    parser.add_argument("--no-barcode", action="store_true", help="Skip barcode extraction (variant calling only)")
    parser.add_argument("--min-base-quality", type=int, default=20, help="Minimum base quality for variant calls (default: 20)")
    parser.add_argument("--expected-barcode-length", type=int, default=15, help="Expected barcode length (default: 15)")
    parser.add_argument("--min-barcode-length", type=int, default=13, help="Minimum acceptable barcode length (default: 13)")
    parser.add_argument("--consensus-threshold", type=float, default=0.9, help="Minimum proportion for barcode consensus (default: 0.9)")
    parser.add_argument("--output-dir", default="variant_calling", help="Output directory (default: variant_calling)")
    parser.add_argument("--sample-name", help="Sample name for output files (default: derived from input BAM)")
    parser.add_argument("--no-index", action="store_true", help="Skip BAM indexing")
    parser.add_argument("--save-pngs", action="store_true", help="Save individual PNG plot files")
    parser.add_argument("--top-variants-nucleosome", type=int, default=10, help="Number of top variants for nucleosome plots (default: 10)")
    parser.add_argument("--expected-variant-types", default="snv,ins,del", help="Variant types for expected coverage (default: snv,ins,del)")
    parser.add_argument("--quiet", action="store_true", help="Suppress progress output")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    
    # Determine if barcode extraction is enabled
    extract_barcode_flag = not args.no_barcode
    
    # Validate barcode arguments if barcode extraction is enabled
    if extract_barcode_flag:
        if args.barcode_start is None or args.barcode_end is None:
            logger.error("--barcode-start and --barcode-end are required unless --no-barcode is specified")
            sys.exit(1)
    
    sample_name = args.sample_name if args.sample_name else Path(args.input_bam).stem
    for suffix in ['.sorted', '.aligned', '.merged', '_sorted', '_aligned', '_merged']:
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
    
    output_dir = Path(args.output_dir)
    report_dir = output_dir / "report"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    output_bam = output_dir / f"{sample_name}.tagged.bam"
    
    try:
        stats, ref_seq = process_bam(
            input_bam=args.input_bam, output_bam=str(output_bam), reference_fasta=args.reference_fasta,
            promoter_chrom=args.promoter_chrom, promoter_start=args.promoter_start, promoter_end=args.promoter_end,
            barcode_start=args.barcode_start, barcode_end=args.barcode_end, min_base_quality=args.min_base_quality,
            expected_barcode_length=args.expected_barcode_length, min_barcode_length=args.min_barcode_length,
            extract_barcode=extract_barcode_flag, quiet=args.quiet)
    except Exception as e:
        logger.error(f"Error processing BAM: {e}")
        sys.exit(1)
    
    if not args.no_index:
        if not args.quiet:
            logger.info("Indexing output BAM...")
        try:
            pysam.index(str(output_bam))
        except Exception as e:
            logger.warning(f"Could not index BAM: {e}")
    
    # Only resolve barcodes if barcode extraction was enabled
    barcode_resolution = {}
    if extract_barcode_flag:
        if not args.quiet:
            logger.info("Resolving barcode-variant mapping...")
        barcode_resolution = resolve_barcode_variants(stats.barcode_variant_map, args.consensus_threshold)
        ambiguous_count = sum(1 for r in barcode_resolution.values() if r.is_ambiguous)
        resolved_count = len(barcode_resolution) - ambiguous_count
        if not args.quiet:
            logger.info(f"Resolved {resolved_count:,} barcodes, {ambiguous_count:,} ambiguous")
    
    variant_types = [v.strip() for v in args.expected_variant_types.split(',')]
    promoter_length = args.promoter_end - args.promoter_start + 1
    coverage = calculate_expected_coverage(stats, promoter_length, ref_seq, variant_types)
    
    if not args.quiet:
        logger.info("Writing reports...")
    
    report_prefix = report_dir / sample_name
    write_stats_report(stats, f"{report_prefix}_stats.tsv")
    write_variants_report(stats, f"{report_prefix}_variants.tsv")
    if extract_barcode_flag:
        write_barcode_variant_map(stats, barcode_resolution, f"{report_prefix}_barcode_variant_map.tsv")
    write_exclusion_report(stats, f"{report_prefix}_exclusion_report.tsv")
    write_expected_coverage_report(coverage, f"{report_prefix}_expected_coverage.tsv")
    
    if not args.quiet:
        logger.info("Generating plots...")
    figures = generate_plots(stats=stats, barcode_resolution=barcode_resolution, coverage=coverage,
        promoter_start=args.promoter_start, promoter_end=args.promoter_end,
        top_variants_nucleosome=args.top_variants_nucleosome,
        output_prefix=str(report_prefix) if args.save_pngs else None, save_individual_pngs=args.save_pngs)
    
    generate_html_report(stats, barcode_resolution, coverage, figures, f"{report_prefix}_report.html", sample_name, 
                         include_barcode=extract_barcode_flag)
    generate_pdf_report(figures, f"{report_prefix}_report.pdf", sample_name)
    
    if not args.quiet:
        print("\n" + "=" * 60)
        if extract_barcode_flag:
            print("MPRA VARIANT-BARCODE TAGGER SUMMARY")
        else:
            print("MPRA VARIANT TAGGER SUMMARY (NO BARCODE)")
        print("=" * 60)
        print(f"Total reads processed:     {stats.total_reads:>12,}")
        if stats.total_reads > 0:
            print(f"Reads passing filters:     {stats.passed_reads:>12,} ({stats.passed_reads/stats.total_reads*100:.1f}%)")
        print(f"  - Wild-type reads:       {stats.wt_reads:>12,}")
        print(f"  - Single variant reads:  {stats.single_variant_reads:>12,}")
        print(f"  - Multi-variant reads:   {stats.multi_variant_reads:>12,}")
        print(f"Unique variants found:     {len(stats.variant_counts):>12,}")
        if extract_barcode_flag:
            ambiguous_count = sum(1 for r in barcode_resolution.values() if r.is_ambiguous)
            resolved_count = len(barcode_resolution) - ambiguous_count
            print(f"Unique barcodes found:     {len(stats.barcode_variant_map):>12,}")
            print(f"  - Resolved barcodes:     {resolved_count:>12,}")
            print(f"  - Ambiguous barcodes:    {ambiguous_count:>12,}")
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