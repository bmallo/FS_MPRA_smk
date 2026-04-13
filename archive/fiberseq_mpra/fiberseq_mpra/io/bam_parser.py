#!/usr/bin/env python3
"""
BAM Parser for Fiberseq MPRA Analysis

This module handles parsing of BAM files containing Fiber-seq footprint data
with custom tags for footprint positions, sizes, nucleosome counts, and variant information.

BAM Tags Used:
    ns: Footprint starts (0-based query coordinates)
    nl: Footprint lengths
    nc: Nucleosome count
    PV: Promoter variant (JSON array of variant IDs, e.g., ["245:A>G"] or "WT")
    VC: Variant count (0=WT, 1=single variant, 2+=multi-variant)
"""

import json
import logging
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Iterator, Set
from collections import defaultdict
import pysam

# Set up logging
logger = logging.getLogger(__name__)


@dataclass
class Footprint:
    """Represents a single footprint on a read."""
    start: int  # 0-based position on the read
    length: int
    
    @property
    def end(self) -> int:
        """Return the end position (exclusive) of the footprint."""
        return self.start + self.length


@dataclass
class ReadFootprints:
    """Represents all footprint information from a single read."""
    read_name: str
    footprints: List[Footprint]
    nucleosome_count: int
    variant_id: Optional[str]  # None for WT, variant ID string for variants
    variant_count: int  # 0 for WT, 1 for single variant, 2+ for multi-variant
    read_length: int
    
    @property
    def is_wt(self) -> bool:
        """Return True if this read is from a wild-type plasmid."""
        return self.variant_count == 0
    
    @property
    def is_single_variant(self) -> bool:
        """Return True if this read has exactly one variant."""
        return self.variant_count == 1


@dataclass
class VariantData:
    """Container for all reads associated with a specific variant."""
    variant_id: str
    reads: List[ReadFootprints] = field(default_factory=list)
    
    @property
    def read_count(self) -> int:
        return len(self.reads)
    
    def add_read(self, read: ReadFootprints):
        self.reads.append(read)


def parse_variant_tag(pv_value: str) -> Tuple[Optional[str], bool]:
    """
    Parse the PV (Promoter Variant) tag value.
    
    Parameters:
    -----------
    pv_value : str
        The PV tag value, either "WT" or a JSON array like '["245:A>G"]'
        
    Returns:
    --------
    Tuple[Optional[str], bool]
        (variant_id or None for WT, is_multi_variant)
    """
    if pv_value == "WT":
        return None, False
    
    try:
        # Try parsing as JSON array
        variants = json.loads(pv_value)
        if isinstance(variants, list):
            if len(variants) == 0:
                return None, False
            elif len(variants) == 1:
                return variants[0], False
            else:
                # Multi-variant - return concatenated ID
                return ";".join(sorted(variants)), True
        else:
            # Single string value
            return str(variants), False
    except json.JSONDecodeError:
        # If not valid JSON, treat as a single variant ID string
        return pv_value, False


def extract_footprints_from_read(read: pysam.AlignedSegment) -> Optional[ReadFootprints]:
    """
    Extract footprint information from a single BAM read.
    
    Parameters:
    -----------
    read : pysam.AlignedSegment
        A single read from a BAM file
        
    Returns:
    --------
    Optional[ReadFootprints]
        Footprint data for the read, or None if required tags are missing
    """
    try:
        # Get required tags
        if not read.has_tag('ns') or not read.has_tag('nl'):
            logger.debug(f"Read {read.query_name} missing ns or nl tags, skipping")
            return None
        
        # Extract footprint starts and lengths
        fp_starts = read.get_tag('ns')
        fp_lengths = read.get_tag('nl')
        
        # Handle different tag formats (could be array or comma-separated string)
        if isinstance(fp_starts, str):
            fp_starts = [int(x) for x in fp_starts.split(',') if x]
        if isinstance(fp_lengths, str):
            fp_lengths = [int(x) for x in fp_lengths.split(',') if x]
        
        # Convert to lists if they're tuples (pysam returns tuples for array tags)
        fp_starts = list(fp_starts) if fp_starts else []
        fp_lengths = list(fp_lengths) if fp_lengths else []
        
        # Validate matching lengths
        if len(fp_starts) != len(fp_lengths):
            logger.warning(f"Read {read.query_name} has mismatched ns/nl lengths: "
                         f"{len(fp_starts)} vs {len(fp_lengths)}")
            return None
        
        # Create Footprint objects
        footprints = [
            Footprint(start=s, length=l) 
            for s, l in zip(fp_starts, fp_lengths)
        ]
        
        # Get nucleosome count (default to 0 if missing)
        nucleosome_count = read.get_tag('nc') if read.has_tag('nc') else 0
        
        # Get variant information
        variant_count = read.get_tag('VC') if read.has_tag('VC') else 0
        
        variant_id = None
        if read.has_tag('PV'):
            pv_value = read.get_tag('PV')
            variant_id, _ = parse_variant_tag(pv_value)
        
        # Get read length
        read_length = read.query_length or read.infer_query_length() or 0
        
        return ReadFootprints(
            read_name=read.query_name,
            footprints=footprints,
            nucleosome_count=nucleosome_count,
            variant_id=variant_id,
            variant_count=variant_count,
            read_length=read_length
        )
        
    except Exception as e:
        logger.warning(f"Error processing read {read.query_name}: {e}")
        return None


def parse_bam_file(
    bam_path: str,
    nucleosome_range: Optional[Tuple[int, int]] = None,
    require_single_variant: bool = True,
    min_read_length: int = 0,
    region: Optional[str] = None
) -> Iterator[ReadFootprints]:
    """
    Parse a BAM file and yield ReadFootprints objects.
    
    Parameters:
    -----------
    bam_path : str
        Path to the BAM file
    nucleosome_range : Optional[Tuple[int, int]]
        If provided, only include reads with nucleosome count in this range (inclusive)
    require_single_variant : bool
        If True, only include WT reads and single-variant reads (exclude multi-variant)
    min_read_length : int
        Minimum read length to include
    region : Optional[str]
        If provided, only parse reads in this region (e.g., "chr1:1000-2000")
        
    Yields:
    -------
    ReadFootprints
        Footprint data for each passing read
    """
    logger.info(f"Opening BAM file: {bam_path}")
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # Determine iterator based on region
        if region:
            read_iter = bam.fetch(region=region)
        else:
            read_iter = bam.fetch()
        
        total_reads = 0
        passed_reads = 0
        filtered_nucleosome = 0
        filtered_multivariant = 0
        filtered_length = 0
        filtered_missing_tags = 0
        
        for read in read_iter:
            total_reads += 1
            
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            
            # Extract footprint data
            read_fp = extract_footprints_from_read(read)
            
            if read_fp is None:
                filtered_missing_tags += 1
                continue
            
            # Apply filters
            
            # Nucleosome count filter
            if nucleosome_range is not None:
                if not (nucleosome_range[0] <= read_fp.nucleosome_count <= nucleosome_range[1]):
                    filtered_nucleosome += 1
                    continue
            
            # Single variant filter
            if require_single_variant and read_fp.variant_count > 1:
                filtered_multivariant += 1
                continue
            
            # Read length filter
            if read_fp.read_length < min_read_length:
                filtered_length += 1
                continue
            
            passed_reads += 1
            yield read_fp
        
        # Log summary statistics
        logger.info(f"BAM parsing complete:")
        logger.info(f"  Total reads processed: {total_reads:,}")
        logger.info(f"  Reads passing filters: {passed_reads:,}")
        logger.info(f"  Filtered (missing tags): {filtered_missing_tags:,}")
        logger.info(f"  Filtered (nucleosome count): {filtered_nucleosome:,}")
        logger.info(f"  Filtered (multi-variant): {filtered_multivariant:,}")
        logger.info(f"  Filtered (read length): {filtered_length:,}")


def load_and_separate_reads(
    bam_path: str,
    nucleosome_range: Optional[Tuple[int, int]] = None,
    require_single_variant: bool = True,
    min_variant_reads: int = 500,
    min_read_length: int = 0,
    region: Optional[str] = None
) -> Tuple[VariantData, Dict[str, VariantData], Dict[str, int]]:
    """
    Load all reads from a BAM file and separate into WT and variant groups.
    
    Parameters:
    -----------
    bam_path : str
        Path to the BAM file
    nucleosome_range : Optional[Tuple[int, int]]
        If provided, only include reads with nucleosome count in this range (inclusive)
    require_single_variant : bool
        If True, only include WT reads and single-variant reads
    min_variant_reads : int
        Minimum number of reads required for a variant to be included in analysis
    min_read_length : int
        Minimum read length to include
    region : Optional[str]
        If provided, only parse reads in this region
        
    Returns:
    --------
    Tuple[VariantData, Dict[str, VariantData], Dict[str, int]]
        (wt_data, variant_dict, excluded_variants)
        - wt_data: VariantData object containing all WT reads
        - variant_dict: Dictionary mapping variant_id to VariantData objects
        - excluded_variants: Dictionary mapping excluded variant_ids to their read counts
    """
    logger.info("Loading and separating reads by variant...")
    
    # Initialize containers
    wt_data = VariantData(variant_id="WT")
    variant_dict: Dict[str, VariantData] = defaultdict(lambda: VariantData(variant_id=""))
    variant_counts: Dict[str, int] = defaultdict(int)
    
    # First pass: count reads per variant
    for read_fp in parse_bam_file(
        bam_path, 
        nucleosome_range=nucleosome_range,
        require_single_variant=require_single_variant,
        min_read_length=min_read_length,
        region=region
    ):
        if read_fp.is_wt:
            wt_data.add_read(read_fp)
        else:
            variant_id = read_fp.variant_id
            variant_counts[variant_id] += 1
            
            # Temporarily store in dict
            if variant_id not in variant_dict:
                variant_dict[variant_id] = VariantData(variant_id=variant_id)
            variant_dict[variant_id].add_read(read_fp)
    
    # Filter variants by minimum read count
    excluded_variants = {}
    passing_variants = {}
    
    for variant_id, var_data in variant_dict.items():
        if var_data.read_count >= min_variant_reads:
            passing_variants[variant_id] = var_data
        else:
            excluded_variants[variant_id] = var_data.read_count
    
    # Log summary
    logger.info(f"Read separation complete:")
    logger.info(f"  WT reads: {wt_data.read_count:,}")
    logger.info(f"  Total variants found: {len(variant_dict):,}")
    logger.info(f"  Variants passing min_reads threshold ({min_variant_reads}): {len(passing_variants):,}")
    logger.info(f"  Variants excluded (low coverage): {len(excluded_variants):,}")
    
    if excluded_variants:
        # Show a few examples of excluded variants
        examples = list(excluded_variants.items())[:5]
        logger.debug(f"  Example excluded variants: {examples}")
    
    return wt_data, passing_variants, excluded_variants


def get_variant_positions(variant_id: str) -> List[Tuple[int, str, str]]:
    """
    Parse a variant ID to extract position and base change information.
    
    Parameters:
    -----------
    variant_id : str
        Variant ID in format "245:A>G" or "245:A>G;301:C>T" for multi-variants
        
    Returns:
    --------
    List[Tuple[int, str, str]]
        List of (position, ref_base, alt_base) tuples
    """
    variants = []
    
    for var in variant_id.split(';'):
        try:
            pos_str, change = var.split(':')
            position = int(pos_str)
            ref_base, alt_base = change.split('>')
            variants.append((position, ref_base.upper(), alt_base.upper()))
        except (ValueError, IndexError) as e:
            logger.warning(f"Could not parse variant ID '{var}': {e}")
            continue
    
    return variants


def get_bam_statistics(bam_path: str) -> Dict:
    """
    Get basic statistics about a BAM file without fully parsing it.
    
    Parameters:
    -----------
    bam_path : str
        Path to the BAM file
        
    Returns:
    --------
    Dict
        Dictionary containing BAM statistics
    """
    stats = {
        'total_reads': 0,
        'mapped_reads': 0,
        'references': [],
        'has_index': False
    }
    
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            stats['references'] = list(bam.references)
            stats['has_index'] = bam.has_index()
            
            # Count reads (this can be slow for large files)
            for read in bam.fetch():
                stats['total_reads'] += 1
                if not read.is_unmapped:
                    stats['mapped_reads'] += 1
                    
    except Exception as e:
        logger.error(f"Error getting BAM statistics: {e}")
        
    return stats


if __name__ == "__main__":
    # Simple test/demo
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    if len(sys.argv) > 1:
        bam_path = sys.argv[1]
        print(f"Testing BAM parser on: {bam_path}")
        
        # Test loading reads
        wt_data, variants, excluded = load_and_separate_reads(
            bam_path,
            min_variant_reads=100,
            require_single_variant=True
        )
        
        print(f"\nWT reads: {wt_data.read_count}")
        print(f"Variants with sufficient coverage: {len(variants)}")
        print(f"Excluded variants: {len(excluded)}")
        
        # Show some variant examples
        for var_id, var_data in list(variants.items())[:3]:
            print(f"\n  Variant: {var_id}")
            print(f"    Reads: {var_data.read_count}")
            positions = get_variant_positions(var_id)
            print(f"    Positions: {positions}")
    else:
        print("Usage: python bam_parser.py <bam_file>")
