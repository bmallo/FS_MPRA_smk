#!/usr/bin/env python3
"""
03_analyze_library.py — Library-Scale Empirical MPRA Fiber-seq Analysis

Tests every variant in a saturation mutagenesis library against WT using
shared null calibrations at representative coverage depths.

Pipeline step 3: tagged BAM → HDF5 + summary TSV + diagnostic PDF

Author: Ben / Stergachis Lab, University of Washington
"""

import argparse
import csv
import json
import logging
import os
import sys
import time
from datetime import datetime

import h5py
import numpy as np

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(iterable, **kwargs):
        return iterable

# Import shared pipeline utilities
sys.path.insert(0, os.path.dirname(__file__))
from common import (
    VERSION, build_analysis_bins, NUC_MIN_LEN_DEFAULT,
    setup_logging, parse_region, safe_hdf5_name, safe_variant_hdf5_name,
    benjamini_hochberg, derive_sample_name, parse_variant_tag,
    parse_variant_id_fields, get_bam_ref_info, resolve_target_chrom,
    parse_bam, group_variants,
    compute_ground_truth_nc,
    run_variant_testing_parallel,
)


# ============================================================================
# HDF5 Output
# ============================================================================

def write_library_hdf5(output_path, parse_stats, ground_truth_nc,
                       null_results_by_depth, all_variant_results,
                       variant_fdr, bin_labels, ref_length, ref_name,
                       promoter_start, promoter_end, analysis_region, args):
    logging.info(f"Writing HDF5: {output_path}")
    a_start, a_end = analysis_region
    n_var = len(all_variant_results)

    with h5py.File(output_path, 'w') as f:
        # ── Metadata ──────────────────────────────────────────────
        meta = f.create_group('metadata')
        meta.attrs['version'] = VERSION
        meta.attrs['date'] = datetime.now().isoformat()
        meta.attrs['reference_name'] = ref_name
        meta.attrs['reference_length'] = ref_length
        meta.attrs['promoter_start'] = promoter_start
        meta.attrs['promoter_end'] = promoter_end
        meta.attrs['analysis_start'] = a_start
        meta.attrs['analysis_end'] = a_end
        meta.attrs['analysis_length'] = a_end - a_start
        meta.attrs['bin_labels'] = [l.encode() for l in bin_labels]
        meta.attrs['n_variants_tested'] = n_var
        meta.attrs['n_null_iterations'] = args.n_null_iterations
        meta.attrs['cluster_threshold_quantile'] = args.cluster_threshold_quantile
        if args.absolute_delta_threshold is not None:
            meta.attrs['absolute_delta_threshold'] = args.absolute_delta_threshold
        meta.attrs['gap_tolerance'] = args.gap_tolerance
        meta.attrs['merge_distance'] = args.merge_distance
        meta.attrs['min_reads'] = args.min_reads
        if args.min_nuc is not None:
            meta.attrs['min_nuc'] = args.min_nuc
        if args.max_nuc is not None:
            meta.attrs['max_nuc'] = args.max_nuc
        meta.attrs['include_multi_variant'] = args.include_multi_variant
        for k, v in parse_stats.items():
            meta.attrs[f'parse_{k}'] = v

        # ── Ground truth NC ───────────────────────────────────────
        if ground_truth_nc is not None:
            gtn = f.create_group('ground_truth_nc')
            gtn.create_dataset('nc_vals', data=ground_truth_nc['nc_vals'])
            gtn.create_dataset('nc_fracs', data=ground_truth_nc['nc_fracs'])
            gtn.attrs['n_reads'] = ground_truth_nc['n_reads']

        # ── WT Occupancy (top-level, one copy) ────────────────────
        # Use the largest-depth null calibration for the most stable estimate
        if null_results_by_depth:
            max_depth = max(null_results_by_depth.keys())
            wt_grp = f.create_group('wt_occupancy')
            wt_grp.attrs['source_null_depth'] = max_depth
            for label in bin_labels:
                wt_grp.create_dataset(
                    safe_hdf5_name(label),
                    data=null_results_by_depth[max_depth]['wt_occ'][label],
                    compression='gzip')

        # ── Null calibrations (summary stats per depth) ───────────
        # Per-variant nulls (Phase 1): no shared per-depth calibration.
        if null_results_by_depth:
            ncg = f.create_group('null_calibration')
            ncg.attrs['depths'] = sorted(list(null_results_by_depth.keys()))
            for depth, nr in null_results_by_depth.items():
                dg = ncg.create_group(str(depth))
                dg.attrs['coverage_level'] = depth
                dg.attrs['n_iterations'] = nr['n_iterations']
                for label in bin_labels:
                    lg = dg.create_group(safe_hdf5_name(label))
                    lg.create_dataset(
                        'pos_null_mean',
                        data=nr['summary'][label]['pos_null_mean'],
                        compression='gzip')
                    lg.create_dataset(
                        'pos_null_std',
                        data=nr['summary'][label]['pos_null_std'],
                        compression='gzip')

        # ── Summary table (flat arrays, fast to load) ─────────────
        sumg = f.create_group('summary')

        # Parse variant IDs into structured fields
        vids = []
        positions = []
        refs = []
        alts = []
        change_types = []
        hdf5_keys = []
        for vr in all_variant_results:
            vid = vr['variant_id']
            vids.append(vid)
            hdf5_keys.append(safe_variant_hdf5_name(vid))
            pos, ref, alt, ct = parse_variant_id_fields(vid)
            positions.append(pos if pos is not None else -1)
            refs.append(ref if ref is not None else '')
            alts.append(alt if alt is not None else '')
            change_types.append(ct if ct is not None else '')

        sumg.create_dataset('variant_ids',
                            data=[v.encode() for v in vids])
        sumg.create_dataset('hdf5_keys',
                            data=[k.encode() for k in hdf5_keys])
        sumg.create_dataset('positions', data=positions)
        sumg.create_dataset('ref_bases',
                            data=[r.encode() for r in refs])
        sumg.create_dataset('alt_bases',
                            data=[a.encode() for a in alts])
        sumg.create_dataset('change_types',
                            data=[c.encode() for c in change_types])
        sumg.create_dataset('n_reads_raw',
                            data=[vr['n_raw'] for vr in all_variant_results])
        sumg.create_dataset('n_reads',
                            data=[vr['n_nc_matched']
                                  for vr in all_variant_results])
        sumg.create_dataset('null_depth_used',
                            data=[vr['null_depth_used']
                                  for vr in all_variant_results])
        sumg.create_dataset('best_cluster_p',
                            data=[vr['best_cluster_p']
                                  for vr in all_variant_results])
        sumg.create_dataset('variant_fdr_q',
                            data=[vr.get('variant_fdr_q', 1.0)
                                  for vr in all_variant_results])
        # NC-shift readout + MDE
        for fld in ('nc_mean_variant', 'nc_mean_wt', 'nc_delta',
                    'nc_wasserstein', 'nc_shift_p', 'nc_shift_q',
                    'mde_median'):
            sumg.create_dataset(
                fld, data=[float(vr.get(fld, float('nan')))
                           for vr in all_variant_results])

        # Per-bin summary columns
        for label in bin_labels:
            sl = safe_hdf5_name(label)
            sumg.create_dataset(f'{sl}_max_abs_delta',
                data=[vr.get(label, {}).get('max_abs_delta', 0.0)
                      for vr in all_variant_results])
            sumg.create_dataset(f'{sl}_n_sig_positions',
                data=[vr.get(label, {}).get('n_sig_positions_fdr10', 0)
                      for vr in all_variant_results])
            sumg.create_dataset(f'{sl}_n_sig_clusters',
                data=[len(vr.get(label, {}).get('significant_clusters', []))
                      for vr in all_variant_results])

        # ── Cluster table (flat, all clusters across all variants) ─
        cluster_rows = []
        for vi, vr in enumerate(all_variant_results):
            vid = vr['variant_id']
            for label in bin_labels:
                lr = vr.get(label, {})
                for ci, cl in enumerate(lr.get('significant_clusters', [])):
                    cluster_rows.append({
                        'variant_idx': vi,
                        'variant_id': vid,
                        'bin_label': label,
                        'cluster_idx': ci,
                        'start': cl['start'],
                        'end': cl['end'],
                        'abs_start': cl['abs_start'],
                        'abs_end': cl['abs_end'],
                        'width': cl['width'],
                        'sum_abs_delta': cl['sum_abs_delta'],
                        'max_abs_delta': cl['max_abs_delta'],
                        'mean_signed_delta': cl.get('mean_signed_delta', 0.0),
                        'direction': cl.get('direction', ''),
                        'sum_p': cl.get('sum_p', 1.0),
                        'peak_position': cl.get('peak_position', 0),
                    })

        if cluster_rows:
            clg = f.create_group('clusters')
            clg.attrs['n_clusters'] = len(cluster_rows)
            clg.create_dataset('variant_idx',
                data=[r['variant_idx'] for r in cluster_rows])
            clg.create_dataset('variant_ids',
                data=[r['variant_id'].encode() for r in cluster_rows])
            clg.create_dataset('bin_labels',
                data=[r['bin_label'].encode() for r in cluster_rows])
            clg.create_dataset('abs_start',
                data=[r['abs_start'] for r in cluster_rows])
            clg.create_dataset('abs_end',
                data=[r['abs_end'] for r in cluster_rows])
            clg.create_dataset('width',
                data=[r['width'] for r in cluster_rows])
            clg.create_dataset('sum_abs_delta',
                data=[r['sum_abs_delta'] for r in cluster_rows])
            clg.create_dataset('max_abs_delta',
                data=[r['max_abs_delta'] for r in cluster_rows])
            clg.create_dataset('mean_signed_delta',
                data=[r['mean_signed_delta'] for r in cluster_rows])
            clg.create_dataset('direction',
                data=[r['direction'].encode() for r in cluster_rows])
            clg.create_dataset('sum_p',
                data=[r['sum_p'] for r in cluster_rows])
            clg.create_dataset('peak_position',
                data=[r['peak_position'] for r in cluster_rows])
            logging.info(f"  {len(cluster_rows)} total significant clusters "
                         f"across {n_var} variants")

        # ── Per-variant detailed arrays ───────────────────────────
        vg = f.create_group('variants')
        for vi, vr in enumerate(all_variant_results):
            vid = vr['variant_id']
            vid_safe = safe_variant_hdf5_name(vid)
            sg = vg.create_group(vid_safe)
            sg.attrs['variant_id'] = vid
            sg.attrs['variant_idx'] = vi
            sg.attrs['n_raw'] = vr['n_raw']
            sg.attrs['n_nc_matched'] = vr['n_nc_matched']
            sg.attrs['null_depth_used'] = vr['null_depth_used']
            sg.attrs['best_cluster_p'] = vr['best_cluster_p']
            sg.attrs['variant_fdr_q'] = vr.get('variant_fdr_q', 1.0)

            # Parsed fields
            pos, ref, alt, ct = parse_variant_id_fields(vid)
            if pos is not None:
                sg.attrs['position'] = pos
            if ref is not None:
                sg.attrs['ref_base'] = ref
            if alt is not None:
                sg.attrs['alt_base'] = alt
            if ct is not None:
                sg.attrs['change_type'] = ct

            for label in bin_labels:
                if label not in vr:
                    continue
                lr = vr[label]
                lg = sg.create_group(safe_hdf5_name(label))
                for key in ['delta_obs', 'variant_occ', 'empirical_p',
                            'q_values', 'z_scores', 'mde']:
                    lg.create_dataset(key, data=lr[key], compression='gzip')
                lg.attrs['n_sig_positions_fdr10'] = lr['n_sig_positions_fdr10']
                lg.attrs['max_abs_delta'] = lr['max_abs_delta']

                # Cluster refs (indices into the flat cluster table)
                sig_cl = lr.get('significant_clusters', [])
                if sig_cl:
                    cgrp = lg.create_group('significant_clusters')
                    for i, cl in enumerate(sig_cl):
                        cg = cgrp.create_group(f'c{i}')
                        for k, v in cl.items():
                            if isinstance(v, (int, float, str, bool)):
                                cg.attrs[k] = v

    logging.info(f"HDF5 written: {output_path}")


# ============================================================================
# Summary TSV
# ============================================================================

def write_summary_tsv(tsv_path, all_variant_results, bin_labels):
    logging.info(f"Writing TSV: {tsv_path}")
    header = ['variant_id', 'n_reads', 'best_cluster_p', 'variant_fdr_q',
              'nc_delta', 'nc_wasserstein', 'nc_shift_p', 'nc_shift_q',
              'mde_median']
    for label in bin_labels:
        short = safe_hdf5_name(label)
        header.extend([f'{short}_max_abs_delta',
                       f'{short}_n_sig_pos_fdr10',
                       f'{short}_n_sig_clusters',
                       f'{short}_top_cluster_start',
                       f'{short}_top_cluster_end',
                       f'{short}_top_cluster_sum_delta',
                       f'{short}_top_cluster_p',
                       f'{short}_top_cluster_direction'])

    with open(tsv_path, 'w', newline='') as fout:
        writer = csv.writer(fout, delimiter='\t')
        writer.writerow(header)
        for vr in all_variant_results:
            row = [vr['variant_id'], vr['n_nc_matched'],
                   f"{vr['best_cluster_p']:.6f}",
                   f"{vr.get('variant_fdr_q', 1.0):.6f}",
                   f"{vr.get('nc_delta', float('nan')):.4f}",
                   f"{vr.get('nc_wasserstein', float('nan')):.4f}",
                   f"{vr.get('nc_shift_p', 1.0):.6f}",
                   f"{vr.get('nc_shift_q', 1.0):.6f}",
                   f"{vr.get('mde_median', float('nan')):.6f}"]
            for label in bin_labels:
                lr = vr.get(label, {})
                row.append(f"{lr.get('max_abs_delta', 0.0):.6f}")
                row.append(str(lr.get('n_sig_positions_fdr10', 0)))
                sc = lr.get('significant_clusters', [])
                row.append(str(len(sc)))
                if sc:
                    top = max(sc, key=lambda c: c['sum_abs_delta'])
                    row.extend([str(top['abs_start']), str(top['abs_end']),
                                f"{top['sum_abs_delta']:.4f}",
                                f"{top.get('max_sum_p', 1.0):.6f}",
                                top.get('direction', '')])
                else:
                    row.extend(['', '', '', '', ''])
            writer.writerow(row)
    logging.info(f"TSV written: {tsv_path}")


# ============================================================================
# Library Diagnostic PDF
# ============================================================================

def generate_library_pdf(pdf_path, all_variant_results, variant_groups,
                         null_results_by_depth, bin_labels,
                         ref_length, promoter_start, promoter_end,
                         analysis_region):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    BIN_COLORS = {
        'sub-TF_10-19bp': '#66c2a5', 'TF_20-40bp': '#fc8d62',
        'PIC_41-80bp': '#8da0cb', 'Nucleosome_81plusbp': '#e78ac3',
    }
    BIN_SHORT = {
        'sub-TF_10-19bp': 'sub-TF', 'TF_20-40bp': 'TF',
        'PIC_41-80bp': 'PIC', 'Nucleosome_81plusbp': 'NUC',
    }

    a_start, a_end = analysis_region
    x = np.arange(a_start, a_end)

    logging.info(f"Generating library PDF: {pdf_path}")

    with PdfPages(pdf_path) as pdf:
        # ---- Page 1: Library overview ----
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Coverage distribution
        ax = axes[0]
        counts = [len(v) for v in variant_groups.values()]
        ax.hist(counts, bins=50, color='#4a86c8', alpha=0.7, edgecolor='white')
        ax.set_xlabel('Reads per Variant')
        ax.set_ylabel('Count')
        ax.set_title(f'Variant Coverage Distribution '
                     f'(n={len(variant_groups):,})')
        ax.axvline(np.median(counts), color='red', linestyle='--',
                    label=f'median={np.median(counts):.0f}')
        ax.legend()

        # Per-variant tested read count (per-variant nulls: no shared
        # depth grid). Falls back to the legacy depth plot if provided.
        ax = axes[1]
        if null_results_by_depth:
            depths = sorted(null_results_by_depth.keys())
            ax.bar(range(len(depths)), depths, color='#fc8d62', alpha=0.7)
            ax.set_xticks(range(len(depths)))
            ax.set_xticklabels(depths, rotation=45)
            ax.set_xlabel('Depth Index')
            ax.set_ylabel('Coverage Level')
            ax.set_title('Null Calibration Depths')
        else:
            ns = [vr.get('n_nc_matched', 0) for vr in all_variant_results]
            if ns:
                ax.hist(ns, bins=40, color='#fc8d62', alpha=0.7,
                        edgecolor='white')
                ax.axvline(np.median(ns), color='red', linestyle='--',
                           label=f'median={np.median(ns):.0f}')
                ax.legend()
            ax.set_xlabel('Reads per Tested Variant (N)')
            ax.set_ylabel('Count')
            ax.set_title('Per-Variant Null Depth (= variant N)')

        fig.suptitle(f'Library Overview: {len(variant_groups):,} variants',
                     fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ---- Page 2: Volcano-style plots per bin ----
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes_flat = axes.flatten()

        for i, label in enumerate(bin_labels):
            if i >= len(axes_flat):
                break
            ax = axes_flat[i]
            max_deltas = []
            neg_log_ps = []
            for vr in all_variant_results:
                lr = vr.get(label, {})
                md = lr.get('max_abs_delta', 0.0)
                # Use signed max delta for direction
                delta = lr.get('delta_obs', np.array([0]))
                peak_idx = np.argmax(np.abs(delta))
                signed_max = float(delta[peak_idx]) if len(delta) > 0 else 0.0
                max_deltas.append(signed_max)
                bp = vr.get('best_cluster_p', 1.0)
                neg_log_ps.append(-np.log10(max(bp, 1e-10)))

            ax.scatter(max_deltas, neg_log_ps, s=8, alpha=0.5,
                       color=BIN_COLORS.get(label, 'gray'))
            ax.axhline(-np.log10(0.05), color='gray', linestyle='--',
                        linewidth=0.5)
            ax.axvline(0, color='gray', linewidth=0.5)
            ax.set_xlabel('Peak Signed Δ')
            ax.set_ylabel('-log10(best cluster p)')
            ax.set_title(BIN_SHORT.get(label, label))

        fig.suptitle('Variant Effect Volcano Plots', fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ---- Page 3: Positional heatmap of significant effects ----
        fig, axes = plt.subplots(len(bin_labels), 1,
                                 figsize=(14, 3 * len(bin_labels)),
                                 sharex=True)
        if len(bin_labels) == 1:
            axes = [axes]

        # Sort variants by position (extract position from variant ID)
        def extract_position(vid):
            try:
                return int(vid.split(':')[0])
            except (ValueError, IndexError):
                return 0

        sorted_results = sorted(all_variant_results,
                                key=lambda v: extract_position(v['variant_id']))

        for i, label in enumerate(bin_labels):
            ax = axes[i]
            # Build heatmap: rows = variants, cols = positions
            n_var = len(sorted_results)
            if n_var == 0:
                continue
            heatmap = np.zeros((n_var, a_end - a_start), dtype=np.float32)
            for vi, vr in enumerate(sorted_results):
                lr = vr.get(label, {})
                delta = lr.get('delta_obs', np.zeros(a_end - a_start))
                heatmap[vi, :] = delta

            vmax = np.percentile(np.abs(heatmap), 99)
            ax.imshow(heatmap, aspect='auto', cmap='RdBu_r',
                      vmin=-vmax, vmax=vmax,
                      extent=[a_start, a_end, n_var, 0])
            ax.set_ylabel(f'{BIN_SHORT.get(label, label)}\nVariants')
            if i == len(bin_labels) - 1:
                ax.set_xlabel('Position')

        fig.suptitle('Positional Effect Heatmap (Δ)', fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ---- Page 4: FDR summary ----
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # P-value distribution
        ax = axes[0]
        best_ps = [vr['best_cluster_p'] for vr in all_variant_results]
        ax.hist(best_ps, bins=50, color='#8da0cb', alpha=0.7,
                edgecolor='white')
        ax.set_xlabel('Best Cluster P-Value')
        ax.set_ylabel('Count')
        ax.set_title('P-Value Distribution')
        ax.axvline(0.05, color='red', linestyle='--',
                    label=f'p=0.05 ({sum(1 for p in best_ps if p < 0.05)} variants)')
        ax.legend()

        # Q-value distribution
        ax = axes[1]
        qs = [vr.get('variant_fdr_q', 1.0) for vr in all_variant_results]
        ax.hist(qs, bins=50, color='#e78ac3', alpha=0.7, edgecolor='white')
        ax.set_xlabel('Variant FDR Q-Value')
        ax.set_ylabel('Count')
        ax.set_title('FDR Q-Value Distribution')
        for thresh in [0.05, 0.10, 0.20]:
            n_sig = sum(1 for q in qs if q < thresh)
            ax.axvline(thresh, color='gray', linestyle='--', linewidth=0.5)
            ax.text(thresh + 0.01, ax.get_ylim()[1] * 0.9,
                    f'FDR<{thresh}: {n_sig}', fontsize=8)

        fig.suptitle('Cross-Variant Multiple Testing', fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ---- Page 5: Summary table for top hits ----
        fig, ax = plt.subplots(figsize=(14, 8))
        ax.axis('off')

        top_hits = sorted(all_variant_results,
                          key=lambda v: v['best_cluster_p'])[:20]
        if top_hits:
            col_labels = ['Variant', 'N reads', 'FDR q',
                          'sub-TF', 'TF', 'PIC', 'NUC']
            table_data = []
            for vr in top_hits:
                row = [vr['variant_id'], str(vr['n_nc_matched']),
                       f"{vr.get('variant_fdr_q', 1.0):.4f}"]
                for label in bin_labels:
                    lr = vr.get(label, {})
                    sc = lr.get('significant_clusters', [])
                    if sc:
                        top_c = max(sc, key=lambda c: c['sum_abs_delta'])
                        row.append(f"{top_c['direction'][0].upper()} "
                                   f"p={top_c.get('sum_p', 1.0):.4f}")
                    else:
                        row.append('—')
                table_data.append(row)

            table = ax.table(cellText=table_data, colLabels=col_labels,
                             loc='center', cellLoc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 1.4)

        ax.set_title('Top 20 Variants by Cluster P-Value', fontsize=14,
                     pad=20)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    logging.info(f"PDF written: {pdf_path}")


# ============================================================================
# Argument Parsing
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=f'Library-Scale Empirical MPRA Analysis (v{VERSION})',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # I/O
    p.add_argument('--bam', required=True, help='Input tagged BAM')
    p.add_argument('-o', '--output-dir', required=True,
                   help='Output directory for HDF5, TSV, and PDF')
    p.add_argument('-n', '--sample-name', default=None,
                   help='Sample name (default: derived from BAM filename)')

    # Regions
    p.add_argument('--target-region', default=None,
                   help='Target region CHROM:START-END or CHROM')
    p.add_argument('--promoter-region', default=None,
                   help='Promoter region START-END (1-based inclusive)')

    # Analysis bins (FiberHMM tracks; all bounds inclusive, configurable)
    p.add_argument('--tf-subtf', type=int, nargs=2, default=[10, 19],
                   metavar=('MIN', 'MAX'),
                   help='tf-track sub_TF length range (default: 10 19)')
    p.add_argument('--tf-tf', type=int, nargs=2, default=[20, 39],
                   metavar=('MIN', 'MAX'),
                   help='tf-track TF length range (default: 20 39)')
    p.add_argument('--tf-pic', type=int, nargs=2, default=[40, 60],
                   metavar=('MIN', 'MAX'),
                   help='tf-track PIC length range (default: 40 60)')
    p.add_argument('--nuc-min-len', type=int, default=NUC_MIN_LEN_DEFAULT,
                   help=f'nuc-track min segment length '
                        f'(default: {NUC_MIN_LEN_DEFAULT})')

    # Filtering
    p.add_argument('--nuc-range', type=int, nargs=2, default=None,
                   metavar=('MIN', 'MAX'),
                   help='Nucleosome count range filter (e.g. 10 25)')
    p.add_argument('--min-reads', type=int, default=50,
                   help='Min reads per variant after filtering (default: 50)')
    p.add_argument('--include-multi-variant', action='store_true',
                   help='Include multi-variant reads (assigned to primary)')
    p.add_argument('--variant-list', default=None,
                   help='File with variant IDs to test (one per line)')

    # Null calibration
    p.add_argument('--n-null-iterations', type=int, default=10000,
                   help='Null subsamples B, common to all variants '
                        '(default: 10000; reported-run value)')
    p.add_argument('--coverage-grid', default=None,
                   help='DEPRECATED, ignored (per-variant nulls)')
    p.add_argument('--random-seed', type=int, default=42)
    # Null stratification (compute lever): reuse one null across
    # variants with similar N and NC distribution.
    # DEFAULT OFF: the WT-vs-WT sweep showed stratification is ~2x
    # anti-conservative (NC-representative mismatch, worse at high N);
    # stratify-OFF is exactly calibrated. Opt-in speed mode only,
    # pending the per-member-reference fix (see redesign plan).
    p.add_argument('--null-stratify', dest='null_stratify',
                   action='store_true', default=False,
                   help='OPT-IN speed mode: reuse nulls across similar '
                        'variants. NOT yet calibrated (~2x anti-'
                        'conservative). Default: off.')
    p.add_argument('--no-null-stratify', dest='null_stratify',
                   action='store_false',
                   help='Independent null per variant (default; calibrated)')
    p.add_argument('--null-strata-n-tol', type=float, default=0.10,
                   help='Max relative N difference within a stratum '
                        '(default: 0.10)')
    p.add_argument('--null-strata-nc-dist', type=float, default=0.30,
                   help='Max NC-distribution Wasserstein-1 within a '
                        'stratum (default: 0.30)')

    # Cluster detection
    p.add_argument('--cluster-threshold-quantile', type=float, default=0.95)
    p.add_argument('--absolute-delta-threshold', type=float, default=None)
    p.add_argument('--gap-tolerance', type=int, default=2)
    p.add_argument('--merge-distance', type=int, default=5)
    p.add_argument('--mde-alpha', type=float, default=0.05,
                   help='Alpha for the minimum-detectable-effect track '
                        '(default: 0.05)')

    # Runtime
    p.add_argument('--threads', type=int, default=None)
    p.add_argument('-v', '--verbose', action='store_true')
    p.add_argument('-q', '--quiet', action='store_true')

    return p.parse_args()


# ============================================================================
# Main
# ============================================================================

def main():
    args = parse_args()
    setup_logging(args.verbose, args.quiet)

    if args.threads is None:
        args.threads = min(os.cpu_count() or 1, 192)

    # Derive sample name
    if args.sample_name is None:
        args.sample_name = derive_sample_name(args.bam)

    # Derive output paths from output-dir + sample-name
    os.makedirs(args.output_dir, exist_ok=True)
    args.output = os.path.join(args.output_dir, f'{args.sample_name}.h5')
    args.tsv = os.path.join(args.output_dir, f'{args.sample_name}_summary.tsv')
    args.pdf = os.path.join(args.output_dir, f'{args.sample_name}_results.pdf')

    # Map --nuc-range to min_nuc/max_nuc for internal use
    args.min_nuc = args.nuc_range[0] if args.nuc_range else None
    args.max_nuc = args.nuc_range[1] if args.nuc_range else None

    logging.info(f"03_analyze_library.py v{VERSION}")
    logging.info(f"BAM: {args.bam}")
    logging.info(f"Sample: {args.sample_name}")
    logging.info(f"Output: {args.output_dir}")
    logging.info(f"Threads: {args.threads}")
    start_time = time.time()

    analysis_bins = build_analysis_bins(
        tf_subTF=tuple(args.tf_subtf), tf_TF=tuple(args.tf_tf),
        tf_PIC=tuple(args.tf_pic), nuc_min_len=args.nuc_min_len)
    bin_labels = [b[0] for b in analysis_bins]
    logging.info(f"Analysis bins: {analysis_bins}")

    target_chrom, target_start, target_end = parse_region(
        args.target_region)

    # ==================================================================
    # Phase A: Parse & Setup
    # ==================================================================
    logging.info("=" * 60)
    logging.info("PHASE A: Parse & Setup")
    logging.info("=" * 60)

    # Determine analysis region BEFORE parsing so coverage matrices
    # are restricted to the region of interest (~15x memory reduction).
    ref_info = get_bam_ref_info(args.bam)
    ref_name, ref_length = resolve_target_chrom(ref_info, target_chrom)

    if args.promoter_region:
        ps = args.promoter_region.split('-')
        prom_s, prom_e = int(ps[0]), int(ps[1]) + 1  # +1: input is inclusive, internal is half-open
    elif target_start is not None:
        prom_s, prom_e = target_start, target_end
    else:
        prom_s, prom_e = 0, ref_length
    logging.info(f"Promoter: {prom_s}-{prom_e - 1} (inclusive, {prom_e - prom_s} bp)")
    analysis_region = (prom_s, prom_e)

    rd, ref_length, ref_name, parse_stats = parse_bam(
        args.bam, analysis_bins, target_chrom=target_chrom,
        analysis_region=analysis_region)

    wt_idx = rd.wt_indices
    logging.info(f"WT reads: {len(wt_idx):,}")

    # Ground truth NC
    ground_truth_nc = compute_ground_truth_nc(
        rd, wt_idx, min_nuc=args.min_nuc, max_nuc=args.max_nuc)
    if ground_truth_nc is None:
        logging.error("No valid NC data. Exiting.")
        sys.exit(1)

    # Group variants
    variant_groups = group_variants(
        rd, include_multi=args.include_multi_variant,
        min_reads=args.min_reads,
        min_nuc=args.min_nuc, max_nuc=args.max_nuc)

    if not variant_groups:
        logging.error("No variants passed filters. Exiting.")
        sys.exit(1)

    # Optional variant list filter
    if args.variant_list:
        with open(args.variant_list) as fin:
            allowed = set(line.strip() for line in fin if line.strip())
        variant_groups = {k: v for k, v in variant_groups.items()
                          if k in allowed}
        logging.info(f"After variant list filter: {len(variant_groups)}")

    # ==================================================================
    # Phase B removed: nulls are now PER-VARIANT (built inside
    # run_variant_testing_parallel at each variant's exact N, WT
    # sampled with replacement, NC-matched to the variant, Option-A
    # NC-reweighted reference). The shared coverage-grid is obsolete.
    if args.coverage_grid:
        logging.warning("--coverage-grid is deprecated and ignored "
                        "(per-variant nulls; see docs/stage3_redesign_plan.md)")
    if args.null_stratify:
        logging.warning(
            "--null-stratify is ON: an OPT-IN speed mode that is NOT yet "
            "calibrated (~2x anti-conservative; WT-vs-WT sweep, worse at "
            "high N). Use the default stratify-OFF for trustworthy FDR.")

    # ==================================================================
    # Phase C: Per-Variant Testing (per-variant nulls)
    # ==================================================================
    logging.info("=" * 60)
    logging.info("PHASE C: Per-Variant Testing")
    logging.info("=" * 60)

    n_total = len(variant_groups)
    all_variant_results = run_variant_testing_parallel(
        rd, wt_idx, variant_groups, analysis_region,
        prom_s, prom_e,
        min_nuc=args.min_nuc, max_nuc=args.max_nuc,
        random_seed=args.random_seed,
        n_null_iterations=args.n_null_iterations,
        cluster_threshold_quantile=args.cluster_threshold_quantile,
        absolute_delta_threshold=args.absolute_delta_threshold,
        gap_tolerance=args.gap_tolerance,
        merge_distance=args.merge_distance,
        n_workers=args.threads,
        stratify=args.null_stratify,
        n_tol=args.null_strata_n_tol,
        nc_dist=args.null_strata_nc_dist,
        mde_alpha=args.mde_alpha)

    n_tested = len(all_variant_results)
    logging.info(f"Tested {n_tested} / {n_total} variants")

    # Cross-variant FDR correction
    if all_variant_results:
        best_ps = np.array([vr['best_cluster_p']
                            for vr in all_variant_results])
        fdr_qs = benjamini_hochberg(best_ps)
        for vr, q in zip(all_variant_results, fdr_qs):
            vr['variant_fdr_q'] = float(q)

        n_sig_005 = int(np.sum(fdr_qs < 0.05))
        n_sig_010 = int(np.sum(fdr_qs < 0.10))
        n_sig_020 = int(np.sum(fdr_qs < 0.20))
        logging.info(f"Significant variants: "
                     f"FDR<0.05={n_sig_005}, "
                     f"FDR<0.10={n_sig_010}, "
                     f"FDR<0.20={n_sig_020}")

        # Independent cross-variant FDR for the NC-shift readout
        nc_ps = np.array([vr.get('nc_shift_p', 1.0)
                          for vr in all_variant_results])
        nc_qs = benjamini_hochberg(nc_ps)
        for vr, q in zip(all_variant_results, nc_qs):
            vr['nc_shift_q'] = float(q)
        logging.info(f"NC-shift significant: "
                     f"FDR<0.10={int(np.sum(nc_qs < 0.10))}")

    # Sort by significance
    all_variant_results.sort(key=lambda v: v['best_cluster_p'])

    # ==================================================================
    # Phase D: Output
    # ==================================================================
    logging.info("=" * 60)
    logging.info("PHASE D: Output")
    logging.info("=" * 60)

    # Per-variant nulls: no shared per-depth null to persist.
    write_library_hdf5(
        args.output, parse_stats, ground_truth_nc,
        None, all_variant_results,
        fdr_qs if all_variant_results else np.array([]),
        bin_labels, ref_length, ref_name,
        prom_s, prom_e, analysis_region, args)

    write_summary_tsv(args.tsv, all_variant_results, bin_labels)

    generate_library_pdf(
        args.pdf, all_variant_results, variant_groups,
        None, bin_labels,
        ref_length, prom_s, prom_e, analysis_region)

    # Summary
    elapsed = time.time() - start_time
    logging.info("")
    logging.info("=" * 70)
    logging.info(f"LIBRARY ANALYSIS COMPLETE (v{VERSION})")
    logging.info("=" * 70)
    logging.info(f"  Variants tested: {n_tested}")
    if all_variant_results:
        logging.info(f"  Significant (FDR<0.10): {n_sig_010}")
        top = all_variant_results[0]
        logging.info(f"  Top hit: {top['variant_id']} "
                     f"(p={top['best_cluster_p']:.6f}, "
                     f"q={top.get('variant_fdr_q', 1.0):.6f})")
    logging.info(f"  Runtime: {elapsed:.1f}s ({elapsed/60:.1f}m)")


if __name__ == '__main__':
    main()
