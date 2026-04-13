"""
Analysis modules for Fiberseq MPRA Analysis

This package contains modules for statistical analysis of footprint data.
"""

from .statistics import (
    TestResult,
    fishers_exact_test,
    test_single_position,
    benjamini_hochberg_correction,
    apply_fdr_correction,
    run_differential_analysis,
    run_all_variant_analyses,
    results_to_dataframe,
    all_results_to_dataframe,
    get_significant_results,
    calculate_log2_fold_change,
    estimate_overdispersion,
)

__all__ = [
    'TestResult',
    'fishers_exact_test',
    'test_single_position',
    'benjamini_hochberg_correction',
    'apply_fdr_correction',
    'run_differential_analysis',
    'run_all_variant_analyses',
    'results_to_dataframe',
    'all_results_to_dataframe',
    'get_significant_results',
    'calculate_log2_fold_change',
    'estimate_overdispersion',
]
