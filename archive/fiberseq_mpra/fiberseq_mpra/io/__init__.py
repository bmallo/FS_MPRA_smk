"""
I/O modules for Fiberseq MPRA Analysis

This package contains modules for reading input files (BAM, FASTA)
and writing output files (TSV, HTML reports).
"""

from .bam_parser import (
    Footprint,
    ReadFootprints,
    VariantData,
    parse_bam_file,
    load_and_separate_reads,
    get_variant_positions,
    get_bam_statistics,
)

from .matrix_builder import (
    SizeBin,
    FootprintMatrix,
    DEFAULT_SIZE_BINS,
    build_footprint_matrix,
    build_matrices_from_variant_data,
    create_contingency_table,
    create_binary_contingency_table,
    parse_size_bins_from_config,
)

__all__ = [
    # BAM parser
    'Footprint',
    'ReadFootprints', 
    'VariantData',
    'parse_bam_file',
    'load_and_separate_reads',
    'get_variant_positions',
    'get_bam_statistics',
    # Matrix builder
    'SizeBin',
    'FootprintMatrix',
    'DEFAULT_SIZE_BINS',
    'build_footprint_matrix',
    'build_matrices_from_variant_data',
    'create_contingency_table',
    'create_binary_contingency_table',
    'parse_size_bins_from_config',
]
