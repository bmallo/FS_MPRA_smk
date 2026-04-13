#!/usr/bin/env python3
"""
Command-Line Interface for Fiberseq MPRA Analysis

This module provides the main entry point for running the analysis pipeline.

Usage:
    python -m fiberseq_mpra.cli.main --bam input.bam --output output_dir/
    
    # Or with configuration file:
    python -m fiberseq_mpra.cli.main --bam input.bam --output output_dir/ --config config.yaml
"""

import argparse
import logging
import sys
import os
from pathlib import Path
from datetime import datetime
from typing import Optional

# Set up logging
def setup_logging(log_level: str = "INFO", log_file: Optional[str] = None):
    """Configure logging for the pipeline."""
    level = getattr(logging, log_level.upper(), logging.INFO)
    
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Fiberseq MPRA Footprint Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python -m fiberseq_mpra.cli.main --bam sample.bam --output results/

    # With nucleosome filtering
    python -m fiberseq_mpra.cli.main --bam sample.bam --output results/ \\
        --nuc-min 10 --nuc-max 25

    # With custom configuration
    python -m fiberseq_mpra.cli.main --bam sample.bam --output results/ \\
        --config config.yaml

    # Generate default config file
    python -m fiberseq_mpra.cli.main --generate-config config.yaml
        """
    )
    
    # Input/output arguments
    parser.add_argument('--bam', '-b', type=str,
                       help='Input BAM file with Fiber-seq footprint data')
    parser.add_argument('--reference', '-r', type=str,
                       help='Reference FASTA file (optional, for sequence annotations)')
    parser.add_argument('--output', '-o', type=str,
                       help='Output directory for results')
    parser.add_argument('--config', '-c', type=str,
                       help='Configuration YAML file (optional)')
    
    # Utility arguments
    parser.add_argument('--generate-config', type=str, metavar='PATH',
                       help='Generate a default configuration file and exit')
    
    # Filtering arguments
    parser.add_argument('--min-reads', type=int,
                       help='Minimum variant reads required for analysis (default: 500)')
    parser.add_argument('--nuc-min', type=int,
                       help='Minimum nucleosome count filter')
    parser.add_argument('--nuc-max', type=int,
                       help='Maximum nucleosome count filter')
    parser.add_argument('--pos-start', type=int,
                       help='Start position for analysis region (0-based)')
    parser.add_argument('--pos-end', type=int,
                       help='End position for analysis region (0-based, inclusive)')
    parser.add_argument('--include-multi-variant', action='store_true',
                       help='Include multi-variant reads (default: single-variant only)')
    
    # Statistical arguments
    parser.add_argument('--fdr', type=float,
                       help='FDR threshold for significance (default: 0.05)')
    parser.add_argument('--min-effect-size', type=float,
                       help='Minimum |log2FC| to report (default: 0)')
    
    # Output arguments
    parser.add_argument('--format', type=str, nargs='+', choices=['html', 'tsv', 'pdf', 'json'],
                       help='Output formats (default: html tsv)')
    parser.add_argument('--per-variant', action='store_true',
                       help='Generate separate output files per variant')
    parser.add_argument('--no-figures', action='store_true',
                       help='Skip static figure generation')
    
    # Other arguments
    parser.add_argument('--log-level', type=str, choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO', help='Logging level (default: INFO)')
    parser.add_argument('--log-file', type=str,
                       help='Log file path (default: stdout only)')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')
    
    return parser.parse_args()


def run_analysis(
    bam_path: str,
    output_dir: str,
    config: 'Config',
    reference_path: Optional[str] = None
) -> dict:
    """
    Run the complete analysis pipeline.
    
    Parameters:
    -----------
    bam_path : str
        Path to input BAM file
    output_dir : str
        Output directory path
    config : Config
        Analysis configuration
    reference_path : str, optional
        Path to reference FASTA
        
    Returns:
    --------
    dict
        Summary statistics and output paths
    """
    from fiberseq_mpra.io.bam_parser import load_and_separate_reads
    from fiberseq_mpra.io.matrix_builder import (
        build_matrices_from_variant_data,
        parse_size_bins_from_config,
    )
    from fiberseq_mpra.analysis.statistics import (
        run_all_variant_analyses,
        all_results_to_dataframe,
    )
    
    logger = logging.getLogger(__name__)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save configuration
    config.save(str(output_path / 'config.yaml'))
    
    results_summary = {
        'start_time': datetime.now().isoformat(),
        'bam_path': bam_path,
        'output_dir': output_dir,
    }
    
    # =========================================================================
    # Step 1: Load and separate reads
    # =========================================================================
    logger.info("=" * 60)
    logger.info("Step 1: Loading BAM file and separating reads")
    logger.info("=" * 60)
    
    wt_data, variant_dict, excluded = load_and_separate_reads(
        bam_path=bam_path,
        nucleosome_range=config.nucleosome_range,
        require_single_variant=config.require_single_variant,
        min_variant_reads=config.min_variant_reads,
        min_read_length=config.min_read_length,
    )
    
    results_summary['wt_read_count'] = wt_data.read_count
    results_summary['variant_count'] = len(variant_dict)
    results_summary['excluded_variant_count'] = len(excluded)
    
    if wt_data.read_count == 0:
        logger.error("No WT reads found! Check BAM file and filters.")
        return results_summary
    
    if len(variant_dict) == 0:
        logger.error("No variants with sufficient coverage found!")
        return results_summary
    
    # =========================================================================
    # Step 2: Determine position range
    # =========================================================================
    logger.info("=" * 60)
    logger.info("Step 2: Determining analysis region")
    logger.info("=" * 60)
    
    if config.position_range:
        position_range = config.position_range
        logger.info(f"Using configured position range: {position_range[0]}-{position_range[1]}")
    else:
        # Infer from read lengths
        all_lengths = [r.read_length for r in wt_data.reads]
        median_length = int(sorted(all_lengths)[len(all_lengths)//2])
        position_range = (0, median_length - 1)
        logger.info(f"Inferred position range from read lengths: {position_range[0]}-{position_range[1]}")
    
    results_summary['position_range'] = position_range
    
    # =========================================================================
    # Step 3: Build footprint matrices
    # =========================================================================
    logger.info("=" * 60)
    logger.info("Step 3: Building footprint matrices")
    logger.info("=" * 60)
    
    size_bins = parse_size_bins_from_config(config.size_bins)
    logger.info(f"Size bins: {[str(b) for b in size_bins]}")
    
    wt_matrix, variant_matrices = build_matrices_from_variant_data(
        wt_data=wt_data,
        variant_dict=variant_dict,
        position_range=position_range,
        size_bins=size_bins,
    )
    
    # =========================================================================
    # Step 4: Run statistical analysis
    # =========================================================================
    logger.info("=" * 60)
    logger.info("Step 4: Running differential analysis")
    logger.info("=" * 60)
    
    all_results = run_all_variant_analyses(wt_matrix, variant_matrices)
    
    # Convert to DataFrame
    results_df = all_results_to_dataframe(all_results)
    
    # Count significant results
    significant_count = results_df['significant'].sum() if len(results_df) > 0 else 0
    logger.info(f"Total tests: {len(results_df):,}")
    logger.info(f"Significant results (FDR < {config.fdr_threshold}): {significant_count:,}")
    
    results_summary['total_tests'] = len(results_df)
    results_summary['significant_tests'] = int(significant_count)
    
    # =========================================================================
    # Step 5: Save results
    # =========================================================================
    logger.info("=" * 60)
    logger.info("Step 5: Saving results")
    logger.info("=" * 60)
    
    output_files = []
    
    # Save TSV
    if 'tsv' in config.output_format:
        tsv_path = output_path / 'all_results.tsv'
        results_df.to_csv(tsv_path, sep='\t', index=False)
        logger.info(f"Saved full results to {tsv_path}")
        output_files.append(str(tsv_path))
        
        # Save significant only
        sig_df = results_df[results_df['significant']]
        if len(sig_df) > 0:
            sig_path = output_path / 'significant_results.tsv'
            sig_df.to_csv(sig_path, sep='\t', index=False)
            logger.info(f"Saved significant results to {sig_path}")
            output_files.append(str(sig_path))
    
    # Save per-variant if requested
    if config.per_variant_output:
        variant_dir = output_path / 'variants'
        variant_dir.mkdir(exist_ok=True)
        
        for var_id in all_results.keys():
            # Clean variant ID for filename
            safe_id = var_id.replace(':', '_').replace('>', '_').replace(';', '__')
            var_df = results_df[results_df['variant_id'] == var_id]
            var_path = variant_dir / f'{safe_id}.tsv'
            var_df.to_csv(var_path, sep='\t', index=False)
        
        logger.info(f"Saved per-variant results to {variant_dir}/")
    
    # Save WT baseline matrix
    wt_baseline_path = output_path / 'wt_baseline.tsv'
    wt_matrix.get_rate_matrix().to_csv(wt_baseline_path, sep='\t')
    logger.info(f"Saved WT baseline matrix to {wt_baseline_path}")
    output_files.append(str(wt_baseline_path))
    
    results_summary['output_files'] = output_files
    results_summary['end_time'] = datetime.now().isoformat()
    
    # =========================================================================
    # Step 6: Generate HTML report (if requested)
    # =========================================================================
    if 'html' in config.output_format:
        logger.info("=" * 60)
        logger.info("Step 6: Generating HTML report")
        logger.info("=" * 60)
        
        try:
            from fiberseq_mpra.visualization.html_report import generate_html_report
            
            html_path = output_path / 'report.html'
            generate_html_report(
                results_df=results_df,
                wt_matrix=wt_matrix,
                variant_matrices=variant_matrices,
                config=config,
                output_path=str(html_path),
            )
            logger.info(f"Generated HTML report: {html_path}")
            output_files.append(str(html_path))
        except ImportError as e:
            logger.warning(f"HTML report generation not available: {e}")
        except Exception as e:
            logger.error(f"Failed to generate HTML report: {e}")
    
    logger.info("=" * 60)
    logger.info("Analysis complete!")
    logger.info("=" * 60)
    
    return results_summary


def main():
    """Main entry point."""
    args = parse_args()
    
    # Set up logging
    setup_logging(args.log_level, args.log_file)
    logger = logging.getLogger(__name__)
    
    # Handle config file generation
    if args.generate_config:
        from fiberseq_mpra.config import create_default_config_file
        create_default_config_file(args.generate_config)
        print(f"Default configuration file created: {args.generate_config}")
        return 0
    
    # Validate required arguments
    if not args.bam:
        logger.error("BAM file is required (use --bam)")
        return 1
    
    if not args.output:
        logger.error("Output directory is required (use --output)")
        return 1
    
    if not os.path.exists(args.bam):
        logger.error(f"BAM file not found: {args.bam}")
        return 1
    
    # Load or create configuration
    from fiberseq_mpra.config import Config, load_config, merge_cli_args
    
    if args.config:
        config = load_config(args.config)
    else:
        config = Config()
    
    # Override with CLI arguments
    cli_overrides = {
        'min_variant_reads': args.min_reads,
        'nuc_min': args.nuc_min,
        'nuc_max': args.nuc_max,
        'pos_start': args.pos_start,
        'pos_end': args.pos_end,
        'fdr_threshold': args.fdr,
        'min_effect_size': args.min_effect_size,
        'output_format': args.format,
    }
    
    # Handle boolean flags
    if args.include_multi_variant:
        cli_overrides['require_single_variant'] = False
    if args.per_variant:
        cli_overrides['per_variant_output'] = True
    if args.no_figures:
        cli_overrides['generate_static_figures'] = False
    
    # Remove None values
    cli_overrides = {k: v for k, v in cli_overrides.items() if v is not None}
    
    if cli_overrides:
        config = merge_cli_args(config, cli_overrides)
    
    # Log configuration
    logger.info("Starting Fiberseq MPRA Analysis")
    logger.info(f"BAM file: {args.bam}")
    logger.info(f"Output directory: {args.output}")
    logger.debug(f"Configuration:\n{config}")
    
    # Run analysis
    try:
        summary = run_analysis(
            bam_path=args.bam,
            output_dir=args.output,
            config=config,
            reference_path=args.reference,
        )
        
        # Print summary
        print("\n" + "=" * 60)
        print("Analysis Summary")
        print("=" * 60)
        print(f"WT reads: {summary.get('wt_read_count', 'N/A'):,}")
        print(f"Variants analyzed: {summary.get('variant_count', 'N/A'):,}")
        print(f"Variants excluded (low coverage): {summary.get('excluded_variant_count', 'N/A'):,}")
        print(f"Total statistical tests: {summary.get('total_tests', 'N/A'):,}")
        print(f"Significant results: {summary.get('significant_tests', 'N/A'):,}")
        print(f"Output directory: {summary.get('output_dir', 'N/A')}")
        print("=" * 60)
        
        return 0
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
