#!/usr/bin/env python3
"""
Matrix Builder for Fiberseq MPRA Analysis

This module constructs position × footprint-size matrices from parsed BAM data.
These matrices represent the footprint landscape and are used for statistical comparisons.

Matrix Structure:
    - Rows: Footprint size bins (e.g., TF: 20-49bp, mid: 50-79bp, nucleosome: 80-200bp)
    - Columns: Positions along the promoter sequence
    - Values: Count of footprints overlapping each position in each size bin
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd

from .bam_parser import ReadFootprints, VariantData, Footprint

# Set up logging
logger = logging.getLogger(__name__)


@dataclass
class SizeBin:
    """Represents a footprint size bin."""
    name: str
    min_size: int  # inclusive
    max_size: int  # inclusive
    
    def contains(self, size: int) -> bool:
        """Check if a footprint size falls within this bin."""
        return self.min_size <= size <= self.max_size
    
    def __str__(self) -> str:
        return f"{self.name} ({self.min_size}-{self.max_size}bp)"


# Default size bins based on biological interpretation
DEFAULT_SIZE_BINS = [
    SizeBin("TF", 20, 49),           # Transcription factor-sized
    SizeBin("mid", 50, 79),          # Medium-sized (large TF or small complex)
    SizeBin("nucleosome", 80, 200),  # Nucleosome-sized
]


@dataclass
class FootprintMatrix:
    """
    Container for a footprint count matrix and associated metadata.
    
    Attributes:
    -----------
    matrix : pd.DataFrame
        DataFrame with size bins as rows and positions as columns
    read_count : int
        Number of reads used to build this matrix
    size_bins : List[SizeBin]
        The size bins used for rows
    position_range : Tuple[int, int]
        The (start, end) positions covered by columns
    variant_id : str
        Identifier for the variant ("WT" for wild-type)
    """
    matrix: pd.DataFrame
    read_count: int
    size_bins: List[SizeBin]
    position_range: Tuple[int, int]
    variant_id: str
    
    @property
    def n_positions(self) -> int:
        """Number of positions (columns) in the matrix."""
        return self.matrix.shape[1]
    
    @property
    def n_bins(self) -> int:
        """Number of size bins (rows) in the matrix."""
        return self.matrix.shape[0]
    
    def get_rate_matrix(self) -> pd.DataFrame:
        """
        Return the matrix normalized by read count (footprint rate per read).
        """
        return self.matrix / self.read_count
    
    def get_counts_at_position(self, position: int) -> pd.Series:
        """Get footprint counts at a specific position across all size bins."""
        if position in self.matrix.columns:
            return self.matrix[position]
        else:
            raise ValueError(f"Position {position} not in matrix (range: {self.position_range})")


def assign_footprint_to_bins(
    footprint: Footprint,
    size_bins: List[SizeBin]
) -> Optional[str]:
    """
    Determine which size bin a footprint belongs to.
    
    Parameters:
    -----------
    footprint : Footprint
        The footprint to classify
    size_bins : List[SizeBin]
        Available size bins
        
    Returns:
    --------
    Optional[str]
        Name of the bin, or None if footprint doesn't fit any bin
    """
    for bin in size_bins:
        if bin.contains(footprint.length):
            return bin.name
    return None


def get_footprint_positions(
    footprint: Footprint,
    position_range: Tuple[int, int]
) -> List[int]:
    """
    Get all positions covered by a footprint within the specified range.
    
    Parameters:
    -----------
    footprint : Footprint
        The footprint
    position_range : Tuple[int, int]
        (start, end) positions to consider (inclusive)
        
    Returns:
    --------
    List[int]
        Positions covered by the footprint within the range
    """
    fp_start = footprint.start
    fp_end = footprint.end  # exclusive
    
    # Clip to position range
    start = max(fp_start, position_range[0])
    end = min(fp_end, position_range[1] + 1)  # +1 because fp_end is exclusive
    
    if start >= end:
        return []
    
    return list(range(start, end))


def build_footprint_matrix(
    reads: List[ReadFootprints],
    position_range: Tuple[int, int],
    size_bins: List[SizeBin] = None,
    variant_id: str = "unknown"
) -> FootprintMatrix:
    """
    Build a footprint count matrix from a list of reads.
    
    Parameters:
    -----------
    reads : List[ReadFootprints]
        List of reads with footprint information
    position_range : Tuple[int, int]
        (start, end) positions to analyze (inclusive, 0-based)
    size_bins : List[SizeBin]
        Footprint size bins to use. Defaults to DEFAULT_SIZE_BINS.
    variant_id : str
        Identifier for this set of reads
        
    Returns:
    --------
    FootprintMatrix
        The constructed matrix with counts
    """
    if size_bins is None:
        size_bins = DEFAULT_SIZE_BINS
    
    # Initialize the count matrix
    positions = list(range(position_range[0], position_range[1] + 1))
    bin_names = [bin.name for bin in size_bins]
    
    # Create matrix filled with zeros
    matrix = pd.DataFrame(
        np.zeros((len(bin_names), len(positions)), dtype=np.int64),
        index=bin_names,
        columns=positions
    )
    
    # Count footprints
    total_footprints = 0
    assigned_footprints = 0
    
    for read in reads:
        for footprint in read.footprints:
            total_footprints += 1
            
            # Determine size bin
            bin_name = assign_footprint_to_bins(footprint, size_bins)
            if bin_name is None:
                continue  # Footprint doesn't fit any bin
            
            # Get positions covered by this footprint
            covered_positions = get_footprint_positions(footprint, position_range)
            if not covered_positions:
                continue  # Footprint doesn't overlap analysis region
            
            assigned_footprints += 1
            
            # Increment counts
            for pos in covered_positions:
                matrix.at[bin_name, pos] += 1
    
    logger.debug(f"Matrix for {variant_id}: {total_footprints} total footprints, "
                f"{assigned_footprints} assigned to bins/positions")
    
    return FootprintMatrix(
        matrix=matrix,
        read_count=len(reads),
        size_bins=size_bins,
        position_range=position_range,
        variant_id=variant_id
    )


def build_matrices_from_variant_data(
    wt_data: VariantData,
    variant_dict: Dict[str, VariantData],
    position_range: Tuple[int, int],
    size_bins: List[SizeBin] = None
) -> Tuple[FootprintMatrix, Dict[str, FootprintMatrix]]:
    """
    Build footprint matrices for WT and all variants.
    
    Parameters:
    -----------
    wt_data : VariantData
        Wild-type read data
    variant_dict : Dict[str, VariantData]
        Dictionary of variant data
    position_range : Tuple[int, int]
        (start, end) positions to analyze
    size_bins : List[SizeBin]
        Footprint size bins to use
        
    Returns:
    --------
    Tuple[FootprintMatrix, Dict[str, FootprintMatrix]]
        (wt_matrix, variant_matrices_dict)
    """
    logger.info(f"Building matrices for positions {position_range[0]}-{position_range[1]}")
    
    # Build WT matrix
    logger.info(f"Building WT matrix ({wt_data.read_count:,} reads)...")
    wt_matrix = build_footprint_matrix(
        reads=wt_data.reads,
        position_range=position_range,
        size_bins=size_bins,
        variant_id="WT"
    )
    
    # Build variant matrices
    variant_matrices = {}
    for i, (var_id, var_data) in enumerate(variant_dict.items()):
        logger.info(f"Building matrix for variant {i+1}/{len(variant_dict)}: "
                   f"{var_id} ({var_data.read_count:,} reads)")
        
        variant_matrices[var_id] = build_footprint_matrix(
            reads=var_data.reads,
            position_range=position_range,
            size_bins=size_bins,
            variant_id=var_id
        )
    
    logger.info(f"Matrix building complete: 1 WT + {len(variant_matrices)} variant matrices")
    
    return wt_matrix, variant_matrices


def create_contingency_table(
    wt_matrix: FootprintMatrix,
    var_matrix: FootprintMatrix,
    size_bin: str,
    position: int
) -> Tuple[np.ndarray, Dict]:
    """
    Create a 2x2 contingency table for Fisher's exact test.
    
    Parameters:
    -----------
    wt_matrix : FootprintMatrix
        Wild-type footprint matrix
    var_matrix : FootprintMatrix
        Variant footprint matrix
    size_bin : str
        Name of the size bin to analyze
    position : int
        Position to analyze
        
    Returns:
    --------
    Tuple[np.ndarray, Dict]
        (contingency_table, metadata)
        
        contingency_table is a 2x2 array:
            [[wt_with_fp, wt_without_fp],
             [var_with_fp, var_without_fp]]
             
        metadata contains additional information about the comparison
    """
    # Get counts
    wt_fp_count = int(wt_matrix.matrix.at[size_bin, position])
    var_fp_count = int(var_matrix.matrix.at[size_bin, position])
    
    wt_total = wt_matrix.read_count
    var_total = var_matrix.read_count
    
    # Note: A read can have multiple footprints at a position, but for the
    # contingency table we use the total footprint count vs total reads.
    # This is an approximation - ideally we'd count reads with/without footprints.
    # For now, we use counts, which is valid when counts << reads.
    
    wt_without_fp = wt_total - wt_fp_count
    var_without_fp = var_total - var_fp_count
    
    # Handle edge case where footprint count exceeds read count
    # (multiple footprints per read at same position)
    wt_without_fp = max(0, wt_without_fp)
    var_without_fp = max(0, var_without_fp)
    
    contingency_table = np.array([
        [wt_fp_count, wt_without_fp],
        [var_fp_count, var_without_fp]
    ])
    
    # Calculate rates for metadata
    wt_rate = wt_fp_count / wt_total if wt_total > 0 else 0
    var_rate = var_fp_count / var_total if var_total > 0 else 0
    
    metadata = {
        'wt_fp_count': wt_fp_count,
        'var_fp_count': var_fp_count,
        'wt_total': wt_total,
        'var_total': var_total,
        'wt_rate': wt_rate,
        'var_rate': var_rate,
        'size_bin': size_bin,
        'position': position
    }
    
    return contingency_table, metadata


def create_binary_contingency_table(
    wt_reads: List[ReadFootprints],
    var_reads: List[ReadFootprints],
    size_bins: List[SizeBin],
    target_bin: str,
    position: int
) -> Tuple[np.ndarray, Dict]:
    """
    Create a 2x2 contingency table based on binary presence/absence of footprints.
    
    This is more statistically proper than using counts, as it avoids
    issues with reads having multiple footprints at the same position.
    
    Parameters:
    -----------
    wt_reads : List[ReadFootprints]
        Wild-type reads
    var_reads : List[ReadFootprints]
        Variant reads
    size_bins : List[SizeBin]
        Size bins to use
    target_bin : str
        Name of the size bin to analyze
    position : int
        Position to analyze
        
    Returns:
    --------
    Tuple[np.ndarray, Dict]
        (contingency_table, metadata)
    """
    def count_reads_with_footprint(reads: List[ReadFootprints]) -> int:
        """Count reads that have at least one footprint of target size at position."""
        count = 0
        for read in reads:
            for fp in read.footprints:
                bin_name = assign_footprint_to_bins(fp, size_bins)
                if bin_name != target_bin:
                    continue
                if fp.start <= position < fp.end:
                    count += 1
                    break  # Only count each read once
        return count
    
    wt_with_fp = count_reads_with_footprint(wt_reads)
    var_with_fp = count_reads_with_footprint(var_reads)
    
    wt_total = len(wt_reads)
    var_total = len(var_reads)
    
    wt_without_fp = wt_total - wt_with_fp
    var_without_fp = var_total - var_with_fp
    
    contingency_table = np.array([
        [wt_with_fp, wt_without_fp],
        [var_with_fp, var_without_fp]
    ])
    
    wt_rate = wt_with_fp / wt_total if wt_total > 0 else 0
    var_rate = var_with_fp / var_total if var_total > 0 else 0
    
    metadata = {
        'wt_with_fp': wt_with_fp,
        'var_with_fp': var_with_fp,
        'wt_total': wt_total,
        'var_total': var_total,
        'wt_rate': wt_rate,
        'var_rate': var_rate,
        'size_bin': target_bin,
        'position': position
    }
    
    return contingency_table, metadata


def parse_size_bins_from_config(config: Dict) -> List[SizeBin]:
    """
    Parse size bin configuration from a dictionary.
    
    Parameters:
    -----------
    config : Dict
        Configuration dictionary with format:
        {
            "TF": [20, 49],
            "mid": [50, 79],
            "nucleosome": [80, 200]
        }
        
    Returns:
    --------
    List[SizeBin]
        List of SizeBin objects
    """
    bins = []
    for name, (min_size, max_size) in config.items():
        bins.append(SizeBin(name=name, min_size=min_size, max_size=max_size))
    
    # Sort by min_size for consistent ordering
    bins.sort(key=lambda b: b.min_size)
    
    return bins


if __name__ == "__main__":
    # Simple test/demo
    import sys
    
    logging.basicConfig(level=logging.DEBUG)
    
    # Create some test data
    print("Testing matrix builder...")
    
    # Create mock footprints
    test_footprints = [
        Footprint(start=10, length=30),   # TF at positions 10-39
        Footprint(start=50, length=60),   # mid at positions 50-109
        Footprint(start=100, length=150), # nucleosome at positions 100-249
    ]
    
    # Create a mock read
    mock_read = ReadFootprints(
        read_name="test_read",
        footprints=test_footprints,
        nucleosome_count=5,
        variant_id=None,
        variant_count=0,
        read_length=300
    )
    
    # Build a matrix
    matrix = build_footprint_matrix(
        reads=[mock_read],
        position_range=(0, 200),
        variant_id="test"
    )
    
    print(f"\nMatrix shape: {matrix.matrix.shape}")
    print(f"Read count: {matrix.read_count}")
    print(f"\nMatrix (first 20 positions):")
    print(matrix.matrix.iloc[:, :20])
    
    print("\n\nSize bins used:")
    for bin in matrix.size_bins:
        print(f"  {bin}")
