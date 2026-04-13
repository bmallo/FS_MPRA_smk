#!/usr/bin/env python3
"""
Statistical Testing for Fiberseq MPRA Analysis

This module provides statistical tests for comparing footprint occupancy
between wild-type and variant sequences.

Primary tests:
    - Fisher's exact test for differential footprint occupancy
    - Benjamini-Hochberg FDR correction for multiple testing

Additional utilities:
    - Effect size calculations (log2 fold change, odds ratio)
    - Overdispersion diagnostics
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from scipy import stats

# Set up logging
logger = logging.getLogger(__name__)


@dataclass
class TestResult:
    """
    Result of a single statistical test comparing WT vs Variant.
    
    Attributes:
    -----------
    position : int
        Genomic position being tested
    size_bin : str
        Footprint size bin being tested
    variant_id : str
        Identifier for the variant
    wt_count : int
        Number of footprints in WT
    var_count : int
        Number of footprints in variant
    wt_total : int
        Total WT reads
    var_total : int
        Total variant reads
    wt_rate : float
        Footprint rate in WT (count/total)
    var_rate : float
        Footprint rate in variant
    pvalue : float
        Raw p-value from Fisher's exact test
    pvalue_adj : float
        FDR-adjusted p-value (set after multiple testing correction)
    log2_fc : float
        Log2 fold change (var_rate / wt_rate)
    odds_ratio : float
        Odds ratio from Fisher's exact test
    direction : str
        'gain' if var > wt, 'loss' if var < wt, 'none' if equal
    """
    position: int
    size_bin: str
    variant_id: str
    wt_count: int
    var_count: int
    wt_total: int
    var_total: int
    wt_rate: float
    var_rate: float
    pvalue: float
    pvalue_adj: float = 1.0
    log2_fc: float = 0.0
    odds_ratio: float = 1.0
    direction: str = 'none'
    
    @property
    def is_significant(self) -> bool:
        """Check if result is significant at FDR < 0.05."""
        return self.pvalue_adj < 0.05
    
    @property
    def neglog10_pvalue(self) -> float:
        """Return -log10(p-value) for visualization."""
        if self.pvalue <= 0:
            return 16  # Cap at -log10(1e-16)
        return -np.log10(max(self.pvalue, 1e-16))
    
    @property
    def neglog10_pvalue_adj(self) -> float:
        """Return -log10(adjusted p-value) for visualization."""
        if self.pvalue_adj <= 0:
            return 16
        return -np.log10(max(self.pvalue_adj, 1e-16))


def calculate_log2_fold_change(
    var_rate: float,
    wt_rate: float,
    pseudocount: float = 1e-6
) -> float:
    """
    Calculate log2 fold change between variant and wild-type rates.
    
    Parameters:
    -----------
    var_rate : float
        Footprint rate in variant
    wt_rate : float
        Footprint rate in wild-type
    pseudocount : float
        Small value to add to avoid log(0)
        
    Returns:
    --------
    float
        Log2(var_rate / wt_rate)
    """
    var_rate_adj = var_rate + pseudocount
    wt_rate_adj = wt_rate + pseudocount
    
    return np.log2(var_rate_adj / wt_rate_adj)


def fishers_exact_test(contingency_table: np.ndarray) -> Tuple[float, float]:
    """
    Perform Fisher's exact test on a 2x2 contingency table.
    
    Parameters:
    -----------
    contingency_table : np.ndarray
        2x2 array: [[wt_with, wt_without], [var_with, var_without]]
        
    Returns:
    --------
    Tuple[float, float]
        (odds_ratio, two_sided_pvalue)
    """
    try:
        odds_ratio, pvalue = stats.fisher_exact(contingency_table)
        return odds_ratio, pvalue
    except Exception as e:
        logger.warning(f"Fisher's exact test failed: {e}")
        return 1.0, 1.0


def test_single_position(
    wt_count: int,
    var_count: int,
    wt_total: int,
    var_total: int,
    position: int,
    size_bin: str,
    variant_id: str
) -> TestResult:
    """
    Perform statistical test at a single position/size bin.
    
    Parameters:
    -----------
    wt_count : int
        Footprint count in WT at this position/bin
    var_count : int
        Footprint count in variant at this position/bin
    wt_total : int
        Total WT reads
    var_total : int
        Total variant reads
    position : int
        Position being tested
    size_bin : str
        Size bin being tested
    variant_id : str
        Variant identifier
        
    Returns:
    --------
    TestResult
        Complete test result
    """
    # Calculate rates
    wt_rate = wt_count / wt_total if wt_total > 0 else 0
    var_rate = var_count / var_total if var_total > 0 else 0
    
    # Build contingency table
    # Using the count-based approach (footprint present vs not at each read position)
    # Note: counts can exceed reads if multiple footprints overlap the position
    wt_without = max(0, wt_total - wt_count)
    var_without = max(0, var_total - var_count)
    
    contingency = np.array([
        [wt_count, wt_without],
        [var_count, var_without]
    ])
    
    # Perform Fisher's exact test
    odds_ratio, pvalue = fishers_exact_test(contingency)
    
    # Calculate log2 fold change
    log2_fc = calculate_log2_fold_change(var_rate, wt_rate)
    
    # Determine direction
    if var_rate > wt_rate:
        direction = 'gain'
    elif var_rate < wt_rate:
        direction = 'loss'
    else:
        direction = 'none'
    
    return TestResult(
        position=position,
        size_bin=size_bin,
        variant_id=variant_id,
        wt_count=wt_count,
        var_count=var_count,
        wt_total=wt_total,
        var_total=var_total,
        wt_rate=wt_rate,
        var_rate=var_rate,
        pvalue=pvalue,
        log2_fc=log2_fc,
        odds_ratio=odds_ratio,
        direction=direction
    )


def benjamini_hochberg_correction(pvalues: List[float]) -> List[float]:
    """
    Apply Benjamini-Hochberg FDR correction to a list of p-values.
    
    Parameters:
    -----------
    pvalues : List[float]
        Raw p-values
        
    Returns:
    --------
    List[float]
        FDR-adjusted p-values (q-values)
    """
    n = len(pvalues)
    if n == 0:
        return []
    
    # Create array with original indices
    pvalue_array = np.array(pvalues)
    
    # Handle NaN values
    nan_mask = np.isnan(pvalue_array)
    valid_pvalues = pvalue_array[~nan_mask]
    valid_indices = np.where(~nan_mask)[0]
    
    if len(valid_pvalues) == 0:
        return [np.nan] * n
    
    # Sort p-values
    sorted_indices = np.argsort(valid_pvalues)
    sorted_pvalues = valid_pvalues[sorted_indices]
    
    # Calculate BH adjusted p-values
    m = len(sorted_pvalues)
    ranks = np.arange(1, m + 1)
    adjusted = sorted_pvalues * m / ranks
    
    # Enforce monotonicity (cumulative minimum from the end)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    
    # Cap at 1.0
    adjusted = np.minimum(adjusted, 1.0)
    
    # Unsort to original order
    unsorted_adjusted = np.zeros(m)
    unsorted_adjusted[sorted_indices] = adjusted
    
    # Reconstruct full array with NaN values
    result = np.full(n, np.nan)
    result[valid_indices] = unsorted_adjusted
    
    return result.tolist()


def apply_fdr_correction(results: List[TestResult]) -> List[TestResult]:
    """
    Apply FDR correction to a list of test results.
    
    Parameters:
    -----------
    results : List[TestResult]
        Test results with raw p-values
        
    Returns:
    --------
    List[TestResult]
        Same results with pvalue_adj field filled in
    """
    if not results:
        return results
    
    # Extract p-values
    pvalues = [r.pvalue for r in results]
    
    # Apply correction
    adjusted = benjamini_hochberg_correction(pvalues)
    
    # Update results
    for result, adj_pvalue in zip(results, adjusted):
        result.pvalue_adj = adj_pvalue
    
    return results


def run_differential_analysis(
    wt_matrix: 'FootprintMatrix',
    var_matrix: 'FootprintMatrix',
    variant_id: str
) -> List[TestResult]:
    """
    Run differential footprint analysis between WT and a variant.
    
    Parameters:
    -----------
    wt_matrix : FootprintMatrix
        Wild-type footprint matrix
    var_matrix : FootprintMatrix
        Variant footprint matrix
    variant_id : str
        Identifier for the variant
        
    Returns:
    --------
    List[TestResult]
        Test results for all positions and size bins
    """
    results = []
    
    # Get dimensions
    size_bins = [bin.name for bin in wt_matrix.size_bins]
    positions = list(wt_matrix.matrix.columns)
    
    wt_total = wt_matrix.read_count
    var_total = var_matrix.read_count
    
    # Test each position and size bin
    for size_bin in size_bins:
        for position in positions:
            wt_count = int(wt_matrix.matrix.at[size_bin, position])
            var_count = int(var_matrix.matrix.at[size_bin, position])
            
            result = test_single_position(
                wt_count=wt_count,
                var_count=var_count,
                wt_total=wt_total,
                var_total=var_total,
                position=position,
                size_bin=size_bin,
                variant_id=variant_id
            )
            results.append(result)
    
    # Apply FDR correction
    results = apply_fdr_correction(results)
    
    logger.info(f"Differential analysis for {variant_id}: "
               f"{len(results)} tests, "
               f"{sum(1 for r in results if r.is_significant)} significant at FDR<0.05")
    
    return results


def run_all_variant_analyses(
    wt_matrix: 'FootprintMatrix',
    variant_matrices: Dict[str, 'FootprintMatrix']
) -> Dict[str, List[TestResult]]:
    """
    Run differential analysis for all variants against WT.
    
    Parameters:
    -----------
    wt_matrix : FootprintMatrix
        Wild-type footprint matrix
    variant_matrices : Dict[str, FootprintMatrix]
        Dictionary of variant matrices
        
    Returns:
    --------
    Dict[str, List[TestResult]]
        Dictionary mapping variant_id to list of test results
    """
    all_results = {}
    
    total_variants = len(variant_matrices)
    for i, (var_id, var_matrix) in enumerate(variant_matrices.items()):
        logger.info(f"Analyzing variant {i+1}/{total_variants}: {var_id}")
        
        results = run_differential_analysis(wt_matrix, var_matrix, var_id)
        all_results[var_id] = results
    
    return all_results


def results_to_dataframe(results: List[TestResult]) -> pd.DataFrame:
    """
    Convert a list of TestResult objects to a pandas DataFrame.
    
    Parameters:
    -----------
    results : List[TestResult]
        List of test results
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with one row per test result
    """
    if not results:
        return pd.DataFrame()
    
    data = []
    for r in results:
        data.append({
            'variant_id': r.variant_id,
            'position': r.position,
            'size_bin': r.size_bin,
            'wt_count': r.wt_count,
            'var_count': r.var_count,
            'wt_total': r.wt_total,
            'var_total': r.var_total,
            'wt_rate': r.wt_rate,
            'var_rate': r.var_rate,
            'log2_fc': r.log2_fc,
            'odds_ratio': r.odds_ratio,
            'pvalue': r.pvalue,
            'pvalue_adj': r.pvalue_adj,
            'neglog10_pvalue': r.neglog10_pvalue,
            'neglog10_pvalue_adj': r.neglog10_pvalue_adj,
            'direction': r.direction,
            'significant': r.is_significant
        })
    
    return pd.DataFrame(data)


def all_results_to_dataframe(
    all_results: Dict[str, List[TestResult]]
) -> pd.DataFrame:
    """
    Convert all variant results to a single DataFrame.
    
    Parameters:
    -----------
    all_results : Dict[str, List[TestResult]]
        Dictionary mapping variant_id to test results
        
    Returns:
    --------
    pd.DataFrame
        Combined DataFrame with all results
    """
    dfs = []
    for var_id, results in all_results.items():
        df = results_to_dataframe(results)
        dfs.append(df)
    
    if not dfs:
        return pd.DataFrame()
    
    return pd.concat(dfs, ignore_index=True)


def get_significant_results(
    results: List[TestResult],
    fdr_threshold: float = 0.05,
    min_log2_fc: float = 0.0
) -> List[TestResult]:
    """
    Filter results to only significant hits.
    
    Parameters:
    -----------
    results : List[TestResult]
        Test results
    fdr_threshold : float
        FDR threshold for significance
    min_log2_fc : float
        Minimum absolute log2 fold change
        
    Returns:
    --------
    List[TestResult]
        Filtered results
    """
    return [
        r for r in results
        if r.pvalue_adj < fdr_threshold and abs(r.log2_fc) >= min_log2_fc
    ]


# ============================================================================
# Overdispersion Diagnostics
# ============================================================================

def estimate_overdispersion(
    wt_reads: List['ReadFootprints'],
    position_range: Tuple[int, int],
    size_bin_name: str,
    size_bins: List['SizeBin'],
    n_subsamples: int = 100,
    subsample_fraction: float = 0.5
) -> Dict:
    """
    Estimate overdispersion by comparing observed variance to binomial expectation.
    
    Parameters:
    -----------
    wt_reads : List[ReadFootprints]
        WT reads to subsample
    position_range : Tuple[int, int]
        Position range to analyze
    size_bin_name : str
        Name of size bin to test
    size_bins : List[SizeBin]
        Size bin definitions
    n_subsamples : int
        Number of subsamples to take
    subsample_fraction : float
        Fraction of reads per subsample
        
    Returns:
    --------
    Dict
        Overdispersion diagnostics including dispersion ratio
    """
    from .matrix_builder import build_footprint_matrix
    
    n_reads = len(wt_reads)
    subsample_size = int(n_reads * subsample_fraction)
    
    if subsample_size < 100:
        logger.warning("Insufficient reads for overdispersion analysis")
        return {'dispersion_ratio': np.nan, 'warning': 'insufficient_reads'}
    
    # Sample multiple times and collect rates at each position
    position_rates = {pos: [] for pos in range(position_range[0], position_range[1] + 1)}
    
    for _ in range(n_subsamples):
        # Random subsample
        indices = np.random.choice(n_reads, size=subsample_size, replace=False)
        subsample = [wt_reads[i] for i in indices]
        
        # Build matrix for subsample
        matrix = build_footprint_matrix(
            reads=subsample,
            position_range=position_range,
            size_bins=size_bins,
            variant_id="subsample"
        )
        
        # Collect rates
        for pos in position_rates.keys():
            if pos in matrix.matrix.columns:
                rate = matrix.matrix.at[size_bin_name, pos] / subsample_size
                position_rates[pos].append(rate)
    
    # Calculate observed variance vs expected (binomial) variance
    dispersion_ratios = []
    
    for pos, rates in position_rates.items():
        if len(rates) < 10:
            continue
        
        mean_rate = np.mean(rates)
        observed_var = np.var(rates)
        
        # Expected variance under binomial: p(1-p)/n
        expected_var = mean_rate * (1 - mean_rate) / subsample_size
        
        if expected_var > 0:
            dispersion_ratios.append(observed_var / expected_var)
    
    if not dispersion_ratios:
        return {'dispersion_ratio': np.nan, 'warning': 'calculation_failed'}
    
    median_dispersion = np.median(dispersion_ratios)
    
    result = {
        'dispersion_ratio': median_dispersion,
        'n_positions_tested': len(dispersion_ratios),
        'overdispersed': median_dispersion > 2.0
    }
    
    if median_dispersion > 2.0:
        logger.warning(f"Significant overdispersion detected (ratio={median_dispersion:.2f}). "
                      f"Consider beta-binomial model.")
    
    return result


if __name__ == "__main__":
    # Simple test/demo
    import numpy as np
    
    logging.basicConfig(level=logging.INFO)
    
    print("Testing statistical functions...")
    
    # Test Fisher's exact
    table = np.array([[10, 90], [25, 75]])
    odds, pval = fishers_exact_test(table)
    print(f"\nFisher's exact test:")
    print(f"  Table: {table.tolist()}")
    print(f"  Odds ratio: {odds:.3f}")
    print(f"  P-value: {pval:.4f}")
    
    # Test BH correction
    pvalues = [0.001, 0.01, 0.02, 0.03, 0.05, 0.1, 0.5]
    adjusted = benjamini_hochberg_correction(pvalues)
    print(f"\nBenjamini-Hochberg correction:")
    for p, q in zip(pvalues, adjusted):
        print(f"  {p:.3f} -> {q:.3f}")
    
    # Test single position test
    result = test_single_position(
        wt_count=100,
        var_count=150,
        wt_total=1000,
        var_total=500,
        position=42,
        size_bin="TF",
        variant_id="test:A>G"
    )
    print(f"\nSingle position test:")
    print(f"  WT rate: {result.wt_rate:.3f}")
    print(f"  Var rate: {result.var_rate:.3f}")
    print(f"  Log2 FC: {result.log2_fc:.3f}")
    print(f"  P-value: {result.pvalue:.4f}")
    print(f"  Direction: {result.direction}")
