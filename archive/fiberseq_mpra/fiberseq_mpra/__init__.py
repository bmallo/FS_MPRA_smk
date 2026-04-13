"""
Fiberseq MPRA Analysis Package

A pipeline for analyzing Fiber-seq footprint data from MPRA experiments
to identify how single nucleotide variants affect chromatin architecture.

Main modules:
    - io: BAM parsing and matrix construction
    - analysis: Statistical testing and multiple testing correction
    - visualization: Heatmaps and interactive HTML reports
    - cli: Command-line interface

Typical usage:
    from fiberseq_mpra import run_analysis
    
    results = run_analysis(
        bam_path="sample.bam",
        reference_path="reference.fa",
        output_dir="output/",
        config=config
    )
"""

__version__ = "0.1.0"
__author__ = "Fiberseq MPRA Analysis Team"

from .config import Config, load_config, DEFAULT_CONFIG

__all__ = [
    'Config',
    'load_config', 
    'DEFAULT_CONFIG',
    '__version__',
]
