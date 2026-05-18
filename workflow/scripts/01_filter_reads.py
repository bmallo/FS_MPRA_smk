#!/usr/bin/env python3
"""
Plasmid Fiber-seq Processor

A processor for plasmid Fiber-seq data that:
- Filters reads by m6A methylation percentage
- Separates reads by chromosome (genomic vs plasmid)
- Applies length and alignment filters to plasmid reads
- Adds nucleosome count (nc) and basepairs-per-nucleosome (bn) tags
- Generates QC plots and a comprehensive filtering report

Author: Ben (Stergachis Lab, University of Washington)
"""

import os
import sys
import argparse
import logging
import datetime
import platform
import subprocess
import shlex
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict

import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm


# =============================================================================
# CONFIGURATION AND CONSTANTS
# =============================================================================

@dataclass
class ProcessingConfig:
    """Configuration for BAM processing."""
    genomic_methylation_range: Tuple[float, float] = (0.1, 1.0)
    plasmid_methylation_range: Tuple[float, float] = (0.1, 1.0)
    length_range: Tuple[int, int] = (50, 50)  # (lower_tolerance, upper_tolerance)
    exact_alignment: bool = True
    skip_methylation: bool = False
    write_genomic: bool = False
    verbose: bool = False


@dataclass
class ProcessingStats:
    """Statistics collected during processing."""
    total_input_reads: int = 0
    total_output_reads: int = 0
    
    # Methylation filtering stats
    methylation_passed_genomic: int = 0
    methylation_failed_genomic: int = 0
    methylation_passed_plasmid: int = 0
    methylation_failed_plasmid: int = 0
    
    # Plasmid-specific filtering stats
    plasmid_failed_length: int = 0
    plasmid_failed_alignment: int = 0
    
    # Chromosome stats
    skipped_chromosomes: int = 0  # chrUn, _random
    unknown_chromosomes: int = 0  # Not in reference
    reads_per_chromosome: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    
    # Nucleosome count distribution (for QC plots)
    nucleosome_counts: Dict[str, Dict[int, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    
    # Methylation distribution (for QC plots) — binned histograms, not raw lists
    # Each value is (counts, bin_edges) from numpy.histogram with 50 bins over [0, 100]
    methylation_histograms: Dict[str, 'numpy.ndarray'] = field(default_factory=dict)
    methylation_sums: Dict[str, float] = field(default_factory=lambda: defaultdict(float))
    methylation_counts: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    
    def get_genomic_total(self) -> int:
        return self.methylation_passed_genomic + self.methylation_failed_genomic
    
    def get_plasmid_total(self) -> int:
        return (self.methylation_passed_plasmid + self.methylation_failed_plasmid + 
                self.plasmid_failed_length + self.plasmid_failed_alignment)


# Genomic chromosome patterns (case-insensitive matching)
GENOMIC_CHROMOSOMES = frozenset([
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
    'chr20', 'chr21', 'chr22', 'chrx', 'chry', 'chrm', 'chrmt', 'chrebv',
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', 
    '15', '16', '17', '18', '19', '20', '21', '22', 'x', 'y', 'm', 'mt', 'ebv'
])

# Cache for chromosome classifications
_chromosome_cache: Dict[str, str] = {}


# =============================================================================
# CHROMOSOME CLASSIFICATION
# =============================================================================

def classify_chromosome(chrom_name: str) -> str:
    """
    Classify a chromosome as 'genomic', 'plasmid', or 'skip'.
    
    Args:
        chrom_name: Name of the chromosome
        
    Returns:
        'genomic' for standard chromosomes, 'plasmid' for others, 'skip' for chrUn/_random
    """
    if chrom_name in _chromosome_cache:
        return _chromosome_cache[chrom_name]
    
    chrom_lower = chrom_name.lower()
    
    # Skip patterns
    if chrom_lower.startswith('chrun') or '_random' in chrom_name.lower():
        result = 'skip'
    elif chrom_lower in GENOMIC_CHROMOSOMES:
        result = 'genomic'
    elif any(chrom_lower.startswith(f'{g}_') for g in GENOMIC_CHROMOSOMES if g.startswith('chr')):
        result = 'genomic'
    else:
        result = 'plasmid'
    
    _chromosome_cache[chrom_name] = result
    return result


def precompute_chromosome_classifications(chromosomes: List[str]) -> Dict[str, str]:
    """Pre-classify all chromosomes and populate the cache."""
    classifications = {}
    for chrom in chromosomes:
        classifications[chrom] = classify_chromosome(chrom)
    return classifications


# =============================================================================
# REFERENCE HANDLING
# =============================================================================

def extract_reference_from_bam_header(bam_path: str) -> Optional[str]:
    """Extract reference FASTA path from BAM header @PG tags."""
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            header = bam_file.header.to_dict()
            
            if 'PG' not in header:
                return None
            
            for pg_entry in header['PG']:
                if 'CL' not in pg_entry:
                    continue
                
                command_line = pg_entry['CL']
                program_name = pg_entry.get('PN', '').lower()
                
                if any(tool in program_name for tool in ['pbmm2', 'minimap2', 'bwa', 'bowtie']):
                    reference = _parse_reference_from_command(command_line)
                    if reference:
                        return reference
            
            return None
            
    except Exception as e:
        logging.error(f"Error reading BAM header: {e}")
        return None


def _parse_reference_from_command(command_line: str) -> Optional[str]:
    """Parse reference FASTA path from an alignment command line."""
    try:
        parts = shlex.split(command_line)
        fasta_extensions = ('.fa', '.fasta', '.fna')
        
        for arg in parts:
            if not arg.startswith('-') and any(arg.endswith(ext) for ext in fasta_extensions):
                return arg
        
        return None
        
    except Exception:
        return None


def parse_fasta_file(fasta_path: str) -> Dict[str, int]:
    """Parse FASTA file and extract chromosome names and lengths."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    reference_dict = {}
    skipped = []
    
    with pysam.FastaFile(fasta_path) as fasta:
        for chrom in fasta.references:
            chrom_type = classify_chromosome(chrom)
            
            if chrom_type == 'skip':
                skipped.append(chrom)
                continue
            
            reference_dict[chrom] = fasta.get_reference_length(chrom)
    
    if skipped:
        logging.info(f"Skipped {len(skipped)} chromosomes (chrUn/_random)")
    
    logging.info(f"Parsed {len(reference_dict)} chromosomes from FASTA")
    return reference_dict


def get_reference_dict(args) -> Tuple[Dict[str, int], str]:
    """Get reference dictionary from arguments or BAM header."""
    if args.fasta_file:
        logging.info(f"Using specified FASTA: {args.fasta_file}")
        return parse_fasta_file(args.fasta_file), args.fasta_file
    
    logging.info("Extracting reference from BAM header...")
    fasta_path = extract_reference_from_bam_header(args.input_bam)
    
    if not fasta_path:
        logging.error("Could not extract reference from BAM header")
        logging.error("Please provide --fasta")
        sys.exit(1)
    
    if not os.path.exists(fasta_path):
        logging.error(f"Reference from BAM header not found: {fasta_path}")
        sys.exit(1)
    
    return parse_fasta_file(fasta_path), fasta_path


# =============================================================================
# FIBER-SEQ DATA PROCESSING
# =============================================================================

def calculate_methylation_percentage(read) -> float:
    """Calculate m6A methylation fraction (0.0-1.0) from modified bases.

    Counts m6A calls from the MM/ML tags (accessible via pysam's
    modified_bases) and divides by A+T base count. This matches the
    original pyft-based calculation.
    """
    seq = read.query_sequence
    if not seq:
        return 0.0

    at_count = seq.count('A') + seq.count('T')
    if at_count == 0:
        return 0.0

    m6a_count = 0
    mod_bases = read.modified_bases
    if mod_bases:
        for (base, strand, mod_type), positions in mod_bases.items():
            if mod_type == 'a':  # m6A modification
                m6a_count += len(positions)

    return m6a_count / at_count


def _nuc_count_from_ma(ma_str) -> Optional[int]:
    """Count nucleosome segments in a FiberHMM MA:Z string.

    MA = "<qlen>;nuc...:s-l,...;msp...:...;tf...:..." — count the
    comma-separated segments in the section whose label starts 'nuc'.
    """
    if not ma_str:
        return None
    for seg in ma_str.split(';')[1:]:
        label, sep, body = seg.partition(':')
        if sep and label.split('+', 1)[0].strip() == 'nuc':
            return 0 if not body else body.count(',') + 1
    return None


def calculate_nucleosome_tags(read) -> Tuple[int, float]:
    """Calculate nucleosome count and basepairs-per-nucleosome.

    Nucleosome-count source priority (the legacy 'as' tag is the MSP /
    accessible track, NOT nucleosomes — using it was a bug):
      1. existing 'nc' tag (FiberHMM source already provides it)
      2. FiberHMM 'MA' string nuc section count
      3. legacy 'ns' tag length (ns/nl = nucleosomes in fibertools)
      4. 0
    """
    nuc_count = None
    if read.has_tag('nc'):
        try:
            nuc_count = int(read.get_tag('nc'))
        except (ValueError, TypeError):
            nuc_count = None
    if nuc_count is None and read.has_tag('MA'):
        nuc_count = _nuc_count_from_ma(read.get_tag('MA'))
    if nuc_count is None:
        try:
            nuc_count = len(read.get_tag('ns'))
        except KeyError:
            nuc_count = 0

    seq_len = read.query_length or 0

    if nuc_count > 0 and seq_len > 0:
        bp_per_nuc = seq_len / nuc_count
    else:
        bp_per_nuc = 0.0

    return nuc_count, bp_per_nuc


# =============================================================================
# FILTERING FUNCTIONS
# =============================================================================

def passes_methylation_filter(methylation_pct: float, chrom_type: str, config: ProcessingConfig) -> bool:
    """Check if a read passes the methylation filter."""
    if config.skip_methylation:
        return True
    
    if chrom_type == 'genomic':
        min_val, max_val = config.genomic_methylation_range
    else:
        min_val, max_val = config.plasmid_methylation_range
    
    return min_val <= methylation_pct <= max_val


def passes_length_filter(read, reference_length: int, config: ProcessingConfig) -> bool:
    """Check if a plasmid read passes the length filter."""
    read_length = read.query_length or 0
    lower_bound = reference_length - config.length_range[0]
    upper_bound = reference_length + config.length_range[1]

    return lower_bound <= read_length <= upper_bound


def passes_alignment_filter(read, reference_length: int) -> bool:
    """Check if a plasmid read spans the entire reference."""
    return read.reference_start == 0 and read.reference_end == reference_length


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_bam_single_pass(
    input_bam: str,
    output_folder: str,
    sample_name: str,
    reference_dict: Dict[str, int],
    config: ProcessingConfig
) -> Tuple[ProcessingStats, Dict[str, str]]:
    """
    Process BAM file in a single pass.

    Writes a single merged plasmid BAM ({sample}.plasmid.bam) and optionally
    a genomic BAM ({sample}.genomic.bam) if config.write_genomic is True.

    Returns (stats, output_bams) where output_bams maps labels to paths,
    e.g. {'plasmid': '/path/to/sample.plasmid.bam', 'genomic': '...'}.
    """
    stats = ProcessingStats()
    output_bams: Dict[str, str] = {}

    # Pre-classify chromosomes
    chrom_classifications = precompute_chromosome_classifications(list(reference_dict.keys()))

    # Count input reads
    logging.info("Counting input reads...")
    stats.total_input_reads = count_reads_fast(input_bam)
    logging.info(f"Input BAM contains {stats.total_input_reads:,} reads")

    # Open template header
    with pysam.AlignmentFile(input_bam, "rb") as template:
        header = template.header

    # Set up output writers
    plasmid_path = os.path.join(output_folder, f"{sample_name}.plasmid.bam")
    plasmid_writer = pysam.AlignmentFile(plasmid_path, "wb", header=header)

    genomic_writer = None
    genomic_path = None
    if config.write_genomic:
        genomic_path = os.path.join(output_folder, f"{sample_name}.genomic.bam")
        genomic_writer = pysam.AlignmentFile(genomic_path, "wb", header=header)

    try:
        pysam_bam = pysam.AlignmentFile(input_bam, "rb")

        for read in tqdm(pysam_bam, total=stats.total_input_reads, desc="Processing reads"):
            chrom = read.reference_name

            # Skip unmapped reads (reference_name is None when --unmapped reads are present)
            if chrom is None:
                stats.unknown_chromosomes += 1
                continue

            # Get classification
            chrom_type = chrom_classifications.get(chrom, classify_chromosome(chrom))

            # Skip chrUn/_random
            if chrom_type == 'skip':
                stats.skipped_chromosomes += 1
                continue

            # Skip unknown chromosomes
            if chrom not in reference_dict:
                stats.unknown_chromosomes += 1
                continue

            # Calculate methylation
            methylation_pct = calculate_methylation_percentage(read)

            # Apply methylation filter
            if not passes_methylation_filter(methylation_pct, chrom_type, config):
                if chrom_type == 'genomic':
                    stats.methylation_failed_genomic += 1
                else:
                    stats.methylation_failed_plasmid += 1
                continue

            # Track methylation pass
            if chrom_type == 'genomic':
                stats.methylation_passed_genomic += 1
            else:
                stats.methylation_passed_plasmid += 1

            # Apply plasmid-specific filters
            if chrom_type == 'plasmid':
                ref_length = reference_dict[chrom]

                if not passes_length_filter(read, ref_length, config):
                    stats.plasmid_failed_length += 1
                    continue

                if config.exact_alignment and not passes_alignment_filter(read, ref_length):
                    stats.plasmid_failed_alignment += 1
                    continue

            # Calculate nucleosome tags
            nuc_count, bp_per_nuc = calculate_nucleosome_tags(read)

            # Add custom tags
            read.set_tag("nc", nuc_count, value_type="i")
            read.set_tag("bn", bp_per_nuc, value_type="f")

            # Write to appropriate output
            if chrom_type == 'plasmid':
                plasmid_writer.write(read)
            elif chrom_type == 'genomic' and genomic_writer is not None:
                genomic_writer.write(read)
            else:
                # Genomic read but --write-genomic not set; skip
                continue

            stats.reads_per_chromosome[chrom] += 1
            stats.total_output_reads += 1

            # Track for QC plots
            stats.nucleosome_counts[chrom][nuc_count] += 1
            # Accumulate methylation into histogram bins (50 bins, 0.0-1.0)
            if chrom not in stats.methylation_histograms:
                stats.methylation_histograms[chrom] = np.zeros(50, dtype=np.int64)
            bin_idx = min(int(methylation_pct * 50), 49)  # 50 bins of width 0.02 over [0,1]
            stats.methylation_histograms[chrom][bin_idx] += 1
            stats.methylation_sums[chrom] += methylation_pct
            stats.methylation_counts[chrom] += 1

        pysam_bam.close()
    finally:
        plasmid_writer.close()
        if genomic_writer is not None:
            genomic_writer.close()

    output_bams['plasmid'] = plasmid_path
    if genomic_path:
        output_bams['genomic'] = genomic_path

    return stats, output_bams


def count_reads_fast(bam_path: str) -> int:
    """Count reads in BAM file using the BAM index statistics."""
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            return sum(s.mapped + s.unmapped
                       for s in bam.get_index_statistics())
    except (ValueError, AttributeError, OSError) as e:
        logging.warning(f"Could not count reads from index: {e}")
        return 0


def index_bam_files(bam_paths: List[str]):
    """Index BAM files."""
    logging.info(f"Indexing {len(bam_paths)} BAM files...")
    for path in bam_paths:
        try:
            subprocess.run(["samtools", "index", path], check=True, capture_output=True)
            logging.debug(f"Indexed: {path}")
        except subprocess.CalledProcessError as e:
            logging.warning(f"Failed to index {path}: {e}")


# =============================================================================
# QC PLOTTING FUNCTIONS
# =============================================================================

def generate_qc_plots(
    stats: ProcessingStats,
    output_folder: str,
    sample_name: str,
    reference_dict: Dict[str, int]
) -> List[str]:
    """Generate QC plots for the processing run."""
    figures_folder = os.path.join(output_folder, 'figures')
    os.makedirs(figures_folder, exist_ok=True)
    
    generated_plots = []
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 150
    
    # 1. Chromosome read count plot
    plot_path = plot_chromosome_read_counts(stats, figures_folder, sample_name)
    if plot_path:
        generated_plots.append(plot_path)
    
    # 2. Nucleosome count distribution plots
    for chrom in stats.nucleosome_counts:
        if stats.nucleosome_counts[chrom]:
            plot_path = plot_nucleosome_distribution(
                stats.nucleosome_counts[chrom], chrom, figures_folder, sample_name
            )
            if plot_path:
                generated_plots.append(plot_path)
    
    # 3. Methylation distribution plots
    plot_path = plot_methylation_distribution(stats, figures_folder, sample_name)
    if plot_path:
        generated_plots.append(plot_path)
    
    # 4. Summary dashboard
    plot_path = plot_summary_dashboard(stats, figures_folder, sample_name, reference_dict)
    if plot_path:
        generated_plots.append(plot_path)
    
    logging.info(f"Generated {len(generated_plots)} QC plots")
    return generated_plots


def plot_chromosome_read_counts(stats: ProcessingStats, figures_folder: str, sample_name: str) -> Optional[str]:
    """Create bar plot of read counts per chromosome."""
    if not stats.reads_per_chromosome:
        return None
    
    genomic_counts = {}
    plasmid_counts = {}
    
    for chrom, count in stats.reads_per_chromosome.items():
        if count > 0:
            chrom_type = classify_chromosome(chrom)
            if chrom_type == 'genomic':
                genomic_counts[chrom] = count
            elif chrom_type == 'plasmid':
                plasmid_counts[chrom] = count
    
    has_genomic = len(genomic_counts) > 0
    has_plasmid = len(plasmid_counts) > 0
    
    if has_genomic and has_plasmid:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    else:
        fig, ax = plt.subplots(figsize=(12, 6))
        ax1 = ax2 = ax
    
    if has_genomic:
        chroms = list(genomic_counts.keys())
        counts = list(genomic_counts.values())
        sorted_pairs = sorted(zip(chroms, counts), key=lambda x: _natural_sort_key(x[0]))
        chroms, counts = zip(*sorted_pairs) if sorted_pairs else ([], [])
        
        ax1.bar(range(len(chroms)), counts, color='steelblue', edgecolor='black', alpha=0.8)
        ax1.set_xticks(range(len(chroms)))
        ax1.set_xticklabels(chroms, rotation=45, ha='right', fontsize=10)
        ax1.set_xlabel('Chromosome', fontsize=12)
        ax1.set_ylabel('Read Count', fontsize=12)
        ax1.set_title(f'Genomic Chromosome Read Counts\n{sample_name}', fontsize=14)
    
    if has_plasmid:
        ax_plasmid = ax2 if has_genomic else ax1
        chroms = list(plasmid_counts.keys())
        counts = list(plasmid_counts.values())
        sorted_pairs = sorted(zip(chroms, counts), key=lambda x: x[1], reverse=True)
        chroms, counts = zip(*sorted_pairs) if sorted_pairs else ([], [])
        
        ax_plasmid.bar(range(len(chroms)), counts, color='coral', edgecolor='black', alpha=0.8)
        ax_plasmid.set_xticks(range(len(chroms)))
        ax_plasmid.set_xticklabels(chroms, rotation=45, ha='right', fontsize=10)
        ax_plasmid.set_xlabel('Plasmid', fontsize=12)
        ax_plasmid.set_ylabel('Read Count', fontsize=12)
        ax_plasmid.set_title(f'Plasmid Read Counts\n{sample_name}', fontsize=14)
    
    plt.tight_layout()
    plot_path = os.path.join(figures_folder, f'{sample_name}_read_counts_by_chromosome.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return plot_path


def plot_nucleosome_distribution(
    nuc_counts: Dict[int, int], chrom: str, figures_folder: str, sample_name: str
) -> Optional[str]:
    """Create histogram of nucleosome count distribution for a chromosome."""
    if not nuc_counts:
        return None
    
    nuc_numbers = sorted(nuc_counts.keys())
    read_counts = [nuc_counts[n] for n in nuc_numbers]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    chrom_type = classify_chromosome(chrom)
    color = 'steelblue' if chrom_type == 'genomic' else 'coral'
    
    ax.bar(nuc_numbers, read_counts, color=color, edgecolor='black', alpha=0.8)
    ax.set_xlabel('Nucleosome Count', fontsize=12)
    ax.set_ylabel('Number of Reads', fontsize=12)
    ax.set_title(f'Nucleosome Count Distribution\n{chrom} - {sample_name}', fontsize=14)
    
    total_reads = sum(read_counts)
    if total_reads > 0:
        weighted_avg = sum(n * c for n, c in nuc_counts.items()) / total_reads
        ax.axvline(x=weighted_avg, color='red', linestyle='--', linewidth=2, label=f'Mean: {weighted_avg:.1f}')
        ax.legend(fontsize=10)
    
    plt.tight_layout()
    plot_path = os.path.join(figures_folder, f'{sample_name}_{chrom}_nucleosome_distribution.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return plot_path


def _aggregate_methylation_histograms(stats, chrom_type_filter):
    """Aggregate pre-binned methylation histograms by chromosome type.

    Returns (counts_array, n_total, mean_val) or (None, 0, 0.0).
    """
    total_hist = np.zeros(50, dtype=np.int64)
    total_sum = 0.0
    total_n = 0
    for chrom, hist in stats.methylation_histograms.items():
        if classify_chromosome(chrom) == chrom_type_filter:
            total_hist += hist
            total_sum += stats.methylation_sums.get(chrom, 0.0)
            total_n += stats.methylation_counts.get(chrom, 0)
    if total_n == 0:
        return None, 0, 0.0
    return total_hist, total_n, total_sum / total_n


def _plot_methylation_hist(ax, hist_counts, n_total, mean_val, color, title):
    """Plot a pre-binned methylation histogram on a matplotlib axis."""
    bin_edges = np.linspace(0, 1.0, 51)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax.bar(bin_centers, hist_counts, width=0.02, color=color,
           edgecolor='black', alpha=0.7)
    ax.set_xlabel('Methylation Percentage', fontsize=12)
    ax.set_ylabel('Number of Reads', fontsize=12)
    ax.set_title(f'{title}\n(n={n_total:,})', fontsize=14)
    ax.axvline(x=mean_val, color='red', linestyle='--', linewidth=2,
               label=f'Mean: {mean_val:.3f}')
    ax.legend(fontsize=10)


def plot_methylation_distribution(stats: ProcessingStats, figures_folder: str, sample_name: str) -> Optional[str]:
    """Create histogram of methylation percentage distribution."""
    gen_hist, gen_n, gen_mean = _aggregate_methylation_histograms(stats, 'genomic')
    pla_hist, pla_n, pla_mean = _aggregate_methylation_histograms(stats, 'plasmid')

    if gen_n == 0 and pla_n == 0:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    if gen_hist is not None:
        _plot_methylation_hist(axes[0], gen_hist, gen_n, gen_mean, 'steelblue',
                               f'Genomic Methylation Distribution\n{sample_name}')
    else:
        axes[0].text(0.5, 0.5, 'No genomic reads', ha='center', va='center', fontsize=14)
        axes[0].set_title('Genomic Methylation Distribution', fontsize=14)

    if pla_hist is not None:
        _plot_methylation_hist(axes[1], pla_hist, pla_n, pla_mean, 'coral',
                               f'Plasmid Methylation Distribution\n{sample_name}')
    else:
        axes[1].text(0.5, 0.5, 'No plasmid reads', ha='center', va='center', fontsize=14)
        axes[1].set_title('Plasmid Methylation Distribution', fontsize=14)

    plt.tight_layout()
    plot_path = os.path.join(figures_folder, f'{sample_name}_methylation_distribution.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()

    return plot_path


def plot_summary_dashboard(
    stats: ProcessingStats, figures_folder: str, sample_name: str, reference_dict: Dict[str, int]
) -> Optional[str]:
    """Create a summary dashboard with multiple QC metrics."""
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. Pie chart: Read disposition
    ax1 = fig.add_subplot(gs[0, 0])
    
    passed = stats.total_output_reads
    failed_meth = stats.methylation_failed_genomic + stats.methylation_failed_plasmid
    failed_length = stats.plasmid_failed_length
    failed_align = stats.plasmid_failed_alignment
    other = stats.skipped_chromosomes + stats.unknown_chromosomes
    
    sizes = [passed, failed_meth, failed_length, failed_align, other]
    labels = ['Passed', 'Failed Meth', 'Failed Length', 'Failed Align', 'Skipped/Unknown']
    colors = ['#2ecc71', '#e74c3c', '#f39c12', '#9b59b6', '#95a5a6']
    
    non_zero = [(s, l, c) for s, l, c in zip(sizes, labels, colors) if s > 0]
    if non_zero:
        sizes, labels, colors = zip(*non_zero)
        ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax1.set_title('Read Disposition', fontsize=12, fontweight='bold')
    
    # 2. Bar chart: Genomic vs Plasmid
    ax2 = fig.add_subplot(gs[0, 1])
    
    genomic_reads = sum(stats.reads_per_chromosome.get(c, 0) 
                       for c in reference_dict if classify_chromosome(c) == 'genomic')
    plasmid_reads = sum(stats.reads_per_chromosome.get(c, 0) 
                       for c in reference_dict if classify_chromosome(c) == 'plasmid')
    
    bars = ax2.bar(['Genomic', 'Plasmid'], [genomic_reads, plasmid_reads], 
                   color=['steelblue', 'coral'], edgecolor='black')
    ax2.set_ylabel('Read Count', fontsize=11)
    ax2.set_title('Reads by Type', fontsize=12, fontweight='bold')
    
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height):,}', ha='center', va='bottom', fontsize=10)
    
    # 3. Stats text box
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.axis('off')
    
    stats_text = f"""Processing Summary
{'='*30}
Total Input Reads:  {stats.total_input_reads:,}
Total Output Reads: {stats.total_output_reads:,}
Pass Rate: {stats.total_output_reads/max(stats.total_input_reads,1)*100:.1f}%

Genomic Reads: {genomic_reads:,}
Plasmid Reads: {plasmid_reads:,}

Filtering Results:
  Meth Failed (Genomic): {stats.methylation_failed_genomic:,}
  Meth Failed (Plasmid): {stats.methylation_failed_plasmid:,}
  Length Failed: {stats.plasmid_failed_length:,}
  Align Failed:  {stats.plasmid_failed_alignment:,}
"""
    ax3.text(0.1, 0.9, stats_text, transform=ax3.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 4. Plasmid read counts
    ax4 = fig.add_subplot(gs[1, :2])
    
    plasmid_chroms = {c: stats.reads_per_chromosome.get(c, 0) 
                      for c in reference_dict if classify_chromosome(c) == 'plasmid'}
    plasmid_chroms = {k: v for k, v in plasmid_chroms.items() if v > 0}
    
    if plasmid_chroms:
        sorted_plasmids = sorted(plasmid_chroms.items(), key=lambda x: x[1], reverse=True)
        names, counts = zip(*sorted_plasmids)
        
        ax4.barh(range(len(names)), counts, color='coral', edgecolor='black', alpha=0.8)
        ax4.set_yticks(range(len(names)))
        ax4.set_yticklabels(names)
        ax4.set_xlabel('Read Count', fontsize=11)
        ax4.set_title('Plasmid Read Counts', fontsize=12, fontweight='bold')
        ax4.invert_yaxis()
    else:
        ax4.text(0.5, 0.5, 'No plasmid reads', ha='center', va='center', fontsize=14)
        ax4.set_title('Plasmid Read Counts', fontsize=12, fontweight='bold')
    
    # 5. Nucleosome distribution for top plasmid
    ax5 = fig.add_subplot(gs[1, 2])
    
    if plasmid_chroms:
        top_plasmid = sorted_plasmids[0][0]
        nuc_counts = stats.nucleosome_counts.get(top_plasmid, {})
        
        if nuc_counts:
            nuc_numbers = sorted(nuc_counts.keys())
            read_counts = [nuc_counts[n] for n in nuc_numbers]
            
            ax5.bar(nuc_numbers, read_counts, color='coral', edgecolor='black', alpha=0.8)
            ax5.set_xlabel('Nucleosome Count', fontsize=11)
            ax5.set_ylabel('Reads', fontsize=11)
            ax5.set_title(f'Nuc Distribution: {top_plasmid}', fontsize=12, fontweight='bold')
    
    # 6. Methylation distributions (from pre-binned histograms)
    ax6 = fig.add_subplot(gs[2, :])

    gen_hist, gen_n, gen_mean = _aggregate_methylation_histograms(stats, 'genomic')
    pla_hist, pla_n, pla_mean = _aggregate_methylation_histograms(stats, 'plasmid')
    bin_edges = np.linspace(0, 1.0, 51)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    if gen_n > 0 or pla_n > 0:
        if gen_hist is not None:
            ax6.bar(bin_centers, gen_hist, width=0.02, alpha=0.6, color='steelblue',
                    label=f'Genomic (n={gen_n:,})', edgecolor='black')
        if pla_hist is not None:
            ax6.bar(bin_centers, pla_hist, width=0.02, alpha=0.6, color='coral',
                    label=f'Plasmid (n={pla_n:,})', edgecolor='black')
        ax6.set_xlabel('Methylation Percentage', fontsize=11)
        ax6.set_ylabel('Number of Reads', fontsize=11)
        ax6.set_title('Methylation Distribution by Chromosome Type', fontsize=12, fontweight='bold')
        ax6.legend(fontsize=10)
    
    plt.suptitle(f'QC Summary Dashboard: {sample_name}', fontsize=16, fontweight='bold', y=1.02)
    
    plot_path = os.path.join(figures_folder, f'{sample_name}_qc_dashboard.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return plot_path


def _natural_sort_key(s: str):
    """Key function for natural sorting of chromosome names."""
    import re
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


# =============================================================================
# REPORT GENERATION
# =============================================================================

def generate_report(
    input_bam: str, output_folder: str, sample_name: str, reference_dict: Dict[str, int],
    fasta_source: str, config: ProcessingConfig, stats: ProcessingStats,
    output_bams: Dict[str, str], command_line: str
) -> str:
    """Generate comprehensive filtering report."""
    report_dir = os.path.join(output_folder, 'report')
    os.makedirs(report_dir, exist_ok=True)
    report_path = os.path.join(report_dir, f"{sample_name}_filtering_report.txt")
    
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    genomic_chroms = [c for c in reference_dict if classify_chromosome(c) == 'genomic']
    plasmid_chroms = [c for c in reference_dict if classify_chromosome(c) == 'plasmid']
    
    genomic_reads = sum(stats.reads_per_chromosome.get(c, 0) for c in genomic_chroms)
    plasmid_reads = sum(stats.reads_per_chromosome.get(c, 0) for c in plasmid_chroms)
    
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PLASMID BAM PROCESSOR - FILTERING REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("RUN INFORMATION\n")
        f.write("-" * 40 + "\n")
        f.write(f"Timestamp:        {timestamp}\n")
        f.write(f"Sample Name:      {sample_name}\n")
        f.write(f"Input BAM:        {os.path.basename(input_bam)}\n")
        f.write(f"Input Path:       {os.path.abspath(input_bam)}\n")
        f.write(f"Output Directory: {os.path.abspath(output_folder)}\n")
        f.write(f"Reference:        {fasta_source}\n")
        f.write(f"Host:             {platform.node()}\n")
        f.write(f"Python:           {platform.python_version()}\n\n")
        
        f.write("FILTERING PARAMETERS\n")
        f.write("-" * 40 + "\n")
        if config.skip_methylation:
            f.write("Methylation Filter: DISABLED\n")
        else:
            f.write(f"Genomic Methylation Range:  {config.genomic_methylation_range[0]:.2f} - {config.genomic_methylation_range[1]:.2f}\n")
            f.write(f"Plasmid Methylation Range:  {config.plasmid_methylation_range[0]:.2f} - {config.plasmid_methylation_range[1]:.2f}\n")
        f.write(f"Plasmid Length Tolerance:   -{config.length_range[0]} / +{config.length_range[1]} bp\n")
        f.write(f"Exact Alignment (Plasmid):  {'Yes' if config.exact_alignment else 'No'}\n")
        f.write(f"Tags Added:                 nc (nucleosome count), bn (bp/nucleosome)\n\n")
        
        f.write("PROCESSING RESULTS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total Input Reads:     {stats.total_input_reads:,}\n")
        f.write(f"Total Output Reads:    {stats.total_output_reads:,}\n")
        if stats.total_input_reads > 0:
            f.write(f"Overall Pass Rate:     {stats.total_output_reads / stats.total_input_reads * 100:.2f}%\n")
        f.write(f"\nGenomic Reads Written: {genomic_reads:,}\n")
        f.write(f"Plasmid Reads Written: {plasmid_reads:,}\n")
        f.write(f"Skipped (chrUn/random):{stats.skipped_chromosomes:,}\n")
        f.write(f"Unknown Chromosomes:   {stats.unknown_chromosomes:,}\n\n")
        
        f.write("FILTERING BREAKDOWN\n")
        f.write("-" * 40 + "\n")
        if not config.skip_methylation:
            f.write("Methylation Filter:\n")
            f.write(f"  Genomic - Passed: {stats.methylation_passed_genomic:,}, Failed: {stats.methylation_failed_genomic:,}\n")
            f.write(f"  Plasmid - Passed: {stats.methylation_passed_plasmid:,}, Failed: {stats.methylation_failed_plasmid:,}\n")
        f.write(f"\nPlasmid-Specific Filters:\n")
        f.write(f"  Failed Length Filter:    {stats.plasmid_failed_length:,}\n")
        f.write(f"  Failed Alignment Filter: {stats.plasmid_failed_alignment:,}\n\n")
        
        f.write("PER-CHROMOSOME BREAKDOWN\n")
        f.write("-" * 40 + "\n")
        f.write(f"{'Chromosome':<20} {'Type':<10} {'Reads':>12} {'Length (bp)':>12}\n")
        f.write("-" * 56 + "\n")
        
        for chrom in sorted(reference_dict.keys()):
            chrom_type = classify_chromosome(chrom)
            if chrom_type == 'skip':
                continue
            reads = stats.reads_per_chromosome.get(chrom, 0)
            length = reference_dict[chrom]
            f.write(f"{chrom:<20} {chrom_type:<10} {reads:>12,} {length:>12,}\n")
        
        f.write("\n")
        
        if output_bams:
            f.write("OUTPUT BAM FILES\n")
            f.write("-" * 40 + "\n")
            for label in sorted(output_bams.keys()):
                rel_path = os.path.relpath(output_bams[label], output_folder)
                f.write(f"{label}: {rel_path}\n")
            f.write("\n")
        
        f.write("COMMAND LINE\n")
        f.write("-" * 40 + "\n")
        f.write(f"{command_line}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("End of Report\n")
        f.write("=" * 80 + "\n")
    
    return report_path


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Process Fiber-seq BAM files: filter by methylation, "
                    "separate genomic/plasmid reads, add nucleosome tags.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (reference from BAM header):
  python 01_filter_reads.py input.bam

  # Specify reference and output location:
  python 01_filter_reads.py input.bam --reference ref.fa -o results/filtered

  # Custom methylation ranges:
  python 01_filter_reads.py input.bam --reference ref.fa \\
      --genomic-methylation-range 0.1 0.8 \\
      --plasmid-methylation-range 0.3 1.0

  # Skip methylation filtering:
  python 01_filter_reads.py input.bam --skip-methylation
"""
    )

    parser.add_argument('input_bam', help='Path to input BAM file')
    parser.add_argument('--reference', dest='fasta_file', default=None,
                        help='Path to reference FASTA (default: extract from BAM header)')

    parser.add_argument('-o', '--output-dir',
                        default=os.path.join(os.getcwd(), 'processed_bams'),
                        help='Output directory (default: ./processed_bams)')
    parser.add_argument('-n', '--sample-name', default=None,
                        help='Sample name (default: derived from BAM filename)')

    parser.add_argument('--genomic-methylation-range', nargs=2, type=float,
                        default=[0.1, 1.0], metavar=('MIN', 'MAX'),
                        help='Methylation range for genomic chromosomes (default: 0.1 1.0)')
    parser.add_argument('--plasmid-methylation-range', nargs=2, type=float,
                        default=[0.1, 1.0], metavar=('MIN', 'MAX'),
                        help='Methylation range for plasmid chromosomes (default: 0.1 1.0)')
    parser.add_argument('--methylation-range', nargs=2, type=float,
                        metavar=('MIN', 'MAX'),
                        help='Methylation range for both types (overrides separate ranges)')
    parser.add_argument('--length-range', nargs=2, type=int, default=[50, 50],
                        metavar=('LOWER', 'UPPER'),
                        help='Length tolerance for plasmid reads (default: 50 50)')
    parser.add_argument('--exact-alignment', action='store_true', default=True,
                        help='Require exact plasmid alignment (default: true)')
    parser.add_argument('--no-exact-alignment', dest='exact_alignment',
                        action='store_false',
                        help='Disable exact alignment requirement')
    parser.add_argument('--skip-methylation', action='store_true',
                        help='Skip methylation filtering')
    parser.add_argument('--write-genomic', action='store_true',
                        help='Also write genomic reads to a separate BAM')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose logging')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress progress output')

    return parser.parse_args()


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def check_dependencies():
    """Verify required tools are available."""
    try:
        subprocess.run(["samtools", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.error("samtools not found. Please install samtools.")
        sys.exit(1)


def main():
    args = parse_arguments()
    setup_logging(args.verbose)

    logging.info("=" * 60)
    logging.info("01_filter_reads — Fiber-seq BAM Processor")
    logging.info("=" * 60)

    check_dependencies()

    if not os.path.exists(args.input_bam):
        logging.error(f"Input BAM not found: {args.input_bam}")
        sys.exit(1)

    # Extract sample name
    if args.sample_name:
        sample_name = args.sample_name
    else:
        sample_name = os.path.basename(args.input_bam)
        if sample_name.endswith('.bam'):
            sample_name = sample_name[:-4]

    logging.info(f"Input: {args.input_bam}")
    logging.info(f"Sample: {sample_name}")
    logging.info(f"Output: {args.output_dir}")

    # Get reference
    reference_dict, fasta_source = get_reference_dict(args)
    logging.info(f"Reference: {fasta_source} ({len(reference_dict)} chromosomes)")

    # Build configuration
    config = ProcessingConfig(
        genomic_methylation_range=tuple(args.methylation_range or args.genomic_methylation_range),
        plasmid_methylation_range=tuple(args.methylation_range or args.plasmid_methylation_range),
        length_range=tuple(args.length_range),
        exact_alignment=args.exact_alignment,
        skip_methylation=args.skip_methylation,
        write_genomic=args.write_genomic,
        verbose=args.verbose
    )

    os.makedirs(args.output_dir, exist_ok=True)

    # Process BAM
    logging.info("Starting processing...")
    stats, output_bams = process_bam_single_pass(
        args.input_bam,
        args.output_dir,
        sample_name,
        reference_dict,
        config
    )

    # Index output BAMs
    if output_bams:
        index_bam_files(list(output_bams.values()))

    # Generate QC plots
    logging.info("Generating QC plots...")
    plot_paths = generate_qc_plots(stats, args.output_dir, sample_name, reference_dict)

    # Generate report
    command_line = ' '.join(sys.argv)
    report_path = generate_report(
        args.input_bam, args.output_dir, sample_name, reference_dict,
        fasta_source, config, stats, output_bams, command_line
    )

    # Print summary
    print("\n" + "=" * 60)
    print("PROCESSING COMPLETE")
    print("=" * 60)
    print(f"Input reads:     {stats.total_input_reads:,}")
    print(f"Output reads:    {stats.total_output_reads:,}")
    if stats.total_input_reads > 0:
        print(f"Pass rate:       {stats.total_output_reads / stats.total_input_reads * 100:.2f}%")
    for label, path in output_bams.items():
        print(f"Output ({label}): {path}")
    print(f"QC plots:        {len(plot_paths)}")
    print(f"Report:          {report_path}")
    print("=" * 60 + "\n")

    logging.info("Done!")
    return 0


if __name__ == '__main__':
    sys.exit(main())