#!/usr/bin/env python3
"""
visualize_calibration.py — Visualization for MPRA Fiber-seq calibration results

Reads HDF5 output from mpra_fiberseq_analysis.py and produces a multi-page
PDF with diagnostic plots organized by analysis phase.

Usage:
    python visualize_calibration.py calibration.h5 [--output plots.pdf]

Dependencies: matplotlib, numpy, h5py
"""

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import seaborn as sns
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.1)
except ImportError:
    pass

# ============================================================================
# Helpers
# ============================================================================

BIN_COLORS = {
    'sub-TF_10-19bp': '#999999',
    'TF_20-40bp': '#e41a1c',
    'PIC_41-80bp': '#377eb8',
    'Nucleosome_81plusbp': '#4daf4a',
}

def get_color(label):
    return BIN_COLORS.get(label, '#333333')

def sn(name):
    """Safe HDF5 name matching the analysis script."""
    return (name.replace(' ', '_').replace('(', '').replace(')', '')
                .replace('+', 'plus').replace('>', 'to').replace(':', '_'))

def load_meta(f):
    m = f['metadata']
    return {
        'version': m.attrs.get('version', '?'),
        'ref_name': m.attrs.get('reference_name', '?'),
        'ref_length': int(m.attrs.get('reference_length', 0)),
        'prom_start': int(m.attrs.get('promoter_start', 0)),
        'prom_end': int(m.attrs.get('promoter_end', 0)),
        'bin_labels': [b.decode() if isinstance(b, bytes) else b
                       for b in m.attrs.get('bin_labels', [])],
        'variant_id': m.attrs.get('variant_id', None),
        'variant_position': m.attrs.get('variant_position', None),
        'parse_total': m.attrs.get('parse_total', 0),
        'parse_parsed': m.attrs.get('parse_parsed', 0),
    }

def vline_variant(ax, vpos, prom_start=0):
    """Draw a vertical line at the variant position."""
    if vpos is not None:
        ax.axvline(vpos - prom_start, color='red', ls='--', lw=1, alpha=0.6, label='Variant')


# ============================================================================
# Phase 1.1: Nucleosome Count Calibration
# ============================================================================

def plot_nuc_calibration(f, meta, pdf):
    if 'nucleosome_calibration' not in f:
        return
    nc = f['nucleosome_calibration']
    bins = meta['bin_labels']
    vpos = meta['variant_position']

    # --- Page 1: Distribution + Cumulative ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Phase 1.1: Nucleosome Count Calibration", fontsize=14, fontweight='bold')

    # Distribution
    wt_d = {int(k): int(v) for k, v in nc['wt_distribution'].attrs.items()}
    var_d = {int(k): int(v) for k, v in nc['var_distribution'].attrs.items()}
    all_nc = sorted(set(list(wt_d.keys()) + list(var_d.keys())))
    x = np.arange(len(all_nc))
    w = 0.35
    ax = axes[0]
    ax.bar(x - w/2, [wt_d.get(n, 0) for n in all_nc], w, label='WT', color='steelblue', alpha=0.8)
    ax.bar(x + w/2, [var_d.get(n, 0) for n in all_nc], w, label='Variant', color='coral', alpha=0.8)
    ax.set_xticks(x); ax.set_xticklabels(all_nc)
    ax.set_xlabel('Nucleosome Count'); ax.set_ylabel('Reads')
    ax.set_title('NC Distribution (WT vs Variant)')
    ks_s = nc.attrs.get('ks_stat', 0); ks_p = nc.attrs.get('ks_pval', 1)
    ax.legend(title=f'KS={ks_s:.3f}, p={ks_p:.2e}')

    # Cumulative: correlation + effect vs min_nc
    cum = nc['cumulative']
    min_ncs = sorted(int(k) for k in cum.keys())
    ax = axes[1]; ax2 = ax.twinx()
    for label in bins:
        s = sn(label)
        corrs, effs, ncs_c, ncs_e = [], [], [], []
        for mnc in min_ncs:
            g = cum[str(mnc)]
            ck = f'corr_{s}'
            if ck in g.attrs:
                corrs.append(float(g.attrs[ck])); ncs_c.append(mnc)
            ek = f'effect_{s}'
            if ek in g.attrs:
                effs.append(float(g.attrs[ek])); ncs_e.append(mnc)
        c = get_color(label)
        ax.plot(ncs_c, corrs, '-o', color=c, ms=4, label=label, alpha=0.8)
        if effs:
            ax2.plot(ncs_e, effs, '--s', color=c, ms=3, alpha=0.5)
    ax.set_xlabel('Min Nucleosome Count'); ax.set_ylabel('Corr with High-NC Ref (solid)')
    ax2.set_ylabel('Effect at Variant (dashed)')
    ax.set_title('Cumulative NC Cutoff')
    ax.legend(fontsize=7, loc='lower left')
    # Annotate read counts
    for mnc in min_ncs[::max(1, len(min_ncs)//6)]:
        g = cum[str(mnc)]
        n = int(g.attrs.get('wt_n', 0))
        ax.annotate(f'{n:,}', (mnc, ax.get_ylim()[0] + 0.01),
                    fontsize=6, ha='center', alpha=0.6, rotation=45)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    pdf.savefig(fig); plt.close(fig)

    # --- Page 2: Per-stratum occupancy heatmaps ---
    if 'strata' not in nc:
        return
    strata = nc['strata']
    strata_ncs = sorted(int(k) for k in strata.keys())
    if not strata_ncs:
        return

    for label in bins:
        s = sn(label)
        profiles = []
        nc_vals = []
        for ncv in strata_ncs:
            key = f'wt_occ_{s}'
            if key in strata[str(ncv)]:
                profiles.append(strata[str(ncv)][key][:])
                nc_vals.append(ncv)
        if not profiles:
            continue

        fig, ax = plt.subplots(figsize=(12, max(3, len(profiles) * 0.4)))
        mat = np.array(profiles)
        im = ax.imshow(mat, aspect='auto', cmap='YlOrRd', interpolation='nearest',
                       extent=[0, meta['ref_length'], len(nc_vals) - 0.5, -0.5])
        ax.set_yticks(range(len(nc_vals)))
        ax.set_yticklabels([f'nc={n}' for n in nc_vals])
        ax.set_xlabel('Position (bp)')
        ax.set_title(f'{label} — WT Occupancy by NC Stratum')
        plt.colorbar(im, ax=ax, label='Occupancy', shrink=0.8)
        if vpos is not None:
            ax.axvline(vpos, color='cyan', ls='--', lw=1.5, alpha=0.8)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Phase 1.2: Quality Calibration
# ============================================================================

def plot_quality_calibration(f, meta, pdf):
    if 'quality_calibration' not in f:
        return
    qc = f['quality_calibration']
    bins = meta['bin_labels']

    thresholds = sorted(int(k) for k in qc.keys())
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Phase 1.2: Footprint Quality Calibration", fontsize=14, fontweight='bold')

    # Split-half correlation vs quality threshold
    ax = axes[0]
    for label in bins:
        s = sn(label)
        corrs = []
        for t in thresholds:
            ck = f'corr_{s}'
            corrs.append(float(qc[str(t)].attrs.get(ck, 0)))
        ax.plot(thresholds, corrs, '-o', color=get_color(label), ms=5, label=label)
    ax.set_xlabel('Min Quality Threshold'); ax.set_ylabel('Split-Half Correlation')
    ax.set_title('Reproducibility vs Quality Filter')
    ax.legend(fontsize=8)
    ax.set_ylim(bottom=0)

    # Footprints retained vs threshold
    ax = axes[1]
    for label in bins:
        s = sn(label)
        nfps = []
        for t in thresholds:
            nk = f'nfp_{s}'
            nfps.append(int(qc[str(t)].attrs.get(nk, 0)))
        ax.plot(thresholds, nfps, '-o', color=get_color(label), ms=5, label=label)
    ax.set_xlabel('Min Quality Threshold'); ax.set_ylabel('Footprints Passing')
    ax.set_title('Data Retention vs Quality Filter')
    ax.legend(fontsize=8)
    ax.set_yscale('log')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Phase 1.3: Power Analysis
# ============================================================================

def plot_power_analysis(f, meta, pdf):
    if 'power_analysis' not in f:
        return
    pa = f['power_analysis']
    bins = meta['bin_labels']

    depths = sorted(int(k) for k in pa.keys())
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("Phase 1.3: Coverage Power Analysis", fontsize=14, fontweight='bold')

    for label in bins:
        s = sn(label)
        det2, det3, mean_z, fpr2, med_z = [], [], [], [], []
        valid_depths = []
        for d in depths:
            if s not in pa[str(d)]:
                continue
            g = pa[str(d)][s]
            valid_depths.append(d)
            det2.append(float(g.attrs.get('detection_rate_z2', 0)))
            det3.append(float(g.attrs.get('detection_rate_z3', 0)))
            mean_z.append(float(g.attrs.get('mean_z', 0)))
            fpr2.append(float(g.attrs.get('fpr_z2', 0)))

        c = get_color(label)
        axes[0].plot(valid_depths, det2, '-o', color=c, ms=4, label=f'{label} |Z|>2')
        axes[0].plot(valid_depths, det3, '--s', color=c, ms=3, alpha=0.6)

        axes[1].plot(valid_depths, mean_z, '-o', color=c, ms=4, label=label)
        axes[2].plot(valid_depths, fpr2, '-o', color=c, ms=4, label=label)

    axes[0].set_xlabel('Read Depth'); axes[0].set_ylabel('Detection Rate')
    axes[0].set_title('Detection Power (solid=Z>2, dashed=Z>3)')
    axes[0].axhline(0.8, ls=':', color='gray', alpha=0.5)
    axes[0].legend(fontsize=7); axes[0].set_xscale('log')

    axes[1].set_xlabel('Read Depth'); axes[1].set_ylabel('Mean Z-score')
    axes[1].set_title('Mean Z at Variant Position')
    axes[1].axhline(0, ls='-', color='gray', alpha=0.3)
    axes[1].legend(fontsize=7); axes[1].set_xscale('log')

    axes[2].set_xlabel('Read Depth'); axes[2].set_ylabel('False Positive Rate')
    axes[2].set_title('FPR at Background Positions (|Z|>2)')
    axes[2].axhline(0.05, ls=':', color='gray', alpha=0.5, label='5%')
    axes[2].legend(fontsize=7); axes[2].set_xscale('log')

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Phase 1.4: WT-vs-WT Null
# ============================================================================

def plot_wt_null(f, meta, pdf):
    if 'wt_null_analysis' not in f:
        return
    wn = f['wt_null_analysis']
    bins = meta['bin_labels']

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Phase 1.4: WT-vs-WT Null Distribution", fontsize=14, fontweight='bold')

    # Summary bar chart
    ax = axes[0]
    labels_present = [l for l in bins if sn(l) in wn]
    x = np.arange(len(labels_present))
    fz2 = [float(wn[sn(l)].attrs.get('frac_z_gt_2', 0)) for l in labels_present]
    fz3 = [float(wn[sn(l)].attrs.get('frac_z_gt_3', 0)) for l in labels_present]
    w = 0.35
    ax.bar(x - w/2, fz2, w, label='|Z|>2', color='steelblue')
    ax.bar(x + w/2, fz3, w, label='|Z|>3', color='coral')
    ax.set_xticks(x); ax.set_xticklabels(labels_present, fontsize=7, rotation=15)
    ax.set_ylabel('Fraction of Positions'); ax.set_title('Null FPR by Bin')
    ax.axhline(0.05, ls=':', color='gray', alpha=0.5, label='5% expected')
    ax.axhline(0.003, ls=':', color='gray', alpha=0.3)
    ax.legend(fontsize=8)

    # Cluster size distribution
    ax = axes[1]
    for l in labels_present:
        if 'max_cluster_sizes' in wn[sn(l)]:
            cs = wn[sn(l)]['max_cluster_sizes'][:]
            if len(cs) > 0:
                ax.hist(cs, bins=max(1, int(max(cs)) - int(min(cs)) + 1),
                        alpha=0.5, color=get_color(l), label=l, density=True)
                p95 = float(wn[sn(l)].attrs.get('cluster_size_95th', 0))
                ax.axvline(p95, color=get_color(l), ls='--', alpha=0.7)
    ax.set_xlabel('Max Cluster Size (bp)'); ax.set_ylabel('Density')
    ax.set_title('Null Cluster Size Distribution (dashed=95th)')
    ax.legend(fontsize=7)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Phase 1.5: Sequence Context
# ============================================================================

def plot_sequence_context(f, meta, pdf):
    if 'sequence_context' not in f:
        return
    sc = f['sequence_context']
    vpos = meta['variant_position']
    ref_len = meta['ref_length']

    fig, axes = plt.subplots(2, 1, figsize=(13, 7))
    fig.suptitle("Phase 1.5: Sequence Context (A/T Content)", fontsize=14, fontweight='bold')

    # Informable positions
    ax = axes[0]
    if 'informable_positions' in sc:
        info = sc['informable_positions'][:].astype(float)
        ax.fill_between(range(len(info)), info, alpha=0.4, color='steelblue', step='mid')
        ax.set_ylabel('Informable (A/T)')
        ax.set_title(f"A/T Positions (total: {sc.attrs.get('at_fraction', 0):.1%} of reference)")
        if vpos is not None:
            ax.axvline(vpos, color='red', ls='--', lw=1.5, alpha=0.7)
        ax.set_xlim(0, ref_len)

    # A/T density sliding windows
    ax = axes[1]
    window_keys = sorted([k for k in sc.keys() if k.startswith('at_density_w')])
    cmap = plt.cm.viridis
    for i, k in enumerate(window_keys):
        w = k.replace('at_density_w', '')
        d = sc[k][:]
        ax.plot(range(len(d)), d, label=f'window={w}bp',
                color=cmap(i / max(1, len(window_keys) - 1)), alpha=0.8)
    ax.set_xlabel('Position (bp)'); ax.set_ylabel('A/T Fraction')
    ax.set_title('Sliding Window A/T Content')
    if vpos is not None:
        ax.axvline(vpos, color='red', ls='--', lw=1.5, alpha=0.7, label='Variant')
    ax.legend(fontsize=8)
    ax.set_xlim(0, ref_len)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Phase 2: WT Landscape & Binding Sites
# ============================================================================

def plot_wt_landscape(f, meta, pdf):
    if 'wt_landscape' not in f:
        return
    wl = f['wt_landscape']
    bins = meta['bin_labels']
    ref_len = meta['ref_length']
    vpos = meta['variant_position']
    ps, pe = meta['prom_start'], meta['prom_end']

    # One page per bin: raw + smoothed + binding sites
    for label in bins:
        s = sn(label)
        rk = f'raw_{s}'
        sk = f'smooth_{s}'
        if rk not in wl or sk not in wl:
            continue

        raw = wl[rk][:]
        smooth = wl[sk][:]

        fig, ax = plt.subplots(figsize=(13, 4))
        positions = np.arange(ref_len)

        ax.fill_between(positions, raw, alpha=0.2, color=get_color(label), label='Raw')
        ax.plot(positions, smooth, color=get_color(label), lw=2, label='Smoothed')

        # Overlay binding sites
        if 'binding_sites' in f and s in f['binding_sites']:
            bs = f['binding_sites'][s]
            for site_key in bs:
                site = bs[site_key]
                ss = int(site.attrs.get('start', 0))
                se = int(site.attrs.get('end', 0))
                pk = float(site.attrs.get('peak_occupancy', 0))
                cv = site.attrs.get('contains_variant', False)
                color = 'red' if cv else 'orange'
                ax.axvspan(ss, se, alpha=0.15, color=color)
                ax.annotate(f'{ss}-{se}', (ss, pk), fontsize=6, alpha=0.7)

        if vpos is not None:
            ax.axvline(vpos, color='red', ls='--', lw=1.5, alpha=0.6, label='Variant')

        ax.axvspan(ps, pe, alpha=0.05, color='blue', label='Promoter')
        ax.set_xlabel('Position (bp)'); ax.set_ylabel('Occupancy')
        ax.set_title(f'Phase 2: WT Landscape — {label} (n={int(wl.attrs.get("n_reads", 0)):,})')
        ax.legend(fontsize=8, loc='upper right')
        ax.set_xlim(0, ref_len)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Phase 3: Position-Level Test Results
# ============================================================================

def plot_position_tests(f, meta, pdf):
    if 'position_tests' not in f:
        return
    pt = f['position_tests']
    bins = meta['bin_labels']
    ref_len = meta['ref_length']
    vpos = meta['variant_position']
    ps, pe = meta['prom_start'], meta['prom_end']

    for vid_key in pt:
        vg = pt[vid_key]

        # Combined Z-score overview: all bins on one page
        fig, axes = plt.subplots(len(bins), 1, figsize=(14, 3 * len(bins)), sharex=True)
        if len(bins) == 1:
            axes = [axes]
        fig.suptitle(f"Phase 3: Position-Level Z-scores — {vid_key}", fontsize=14, fontweight='bold')

        for i, label in enumerate(bins):
            s = sn(label)
            if s not in vg:
                continue
            lg = vg[s]
            ax = axes[i]

            z = lg['z_scores'][:]
            positions = np.arange(len(z))
            vo = lg['variant_occupancy'][:]
            nm = lg['null_mean'][:]

            # Z-score bar plot
            colors = np.where(z > 0, '#e41a1c', '#377eb8')
            ax.bar(positions, z, width=1.0, color=colors, alpha=0.7, linewidth=0)
            ax.axhline(2, ls=':', color='gray', alpha=0.4)
            ax.axhline(-2, ls=':', color='gray', alpha=0.4)
            ax.axhline(3, ls=':', color='orange', alpha=0.4)
            ax.axhline(-3, ls=':', color='orange', alpha=0.4)

            if vpos is not None:
                ax.axvline(vpos, color='red', ls='--', lw=2, alpha=0.7)

            # Mark significant clusters
            for zt in [2, 3]:
                ck = f'clusters_z{zt}'
                if ck in lg:
                    for c_key in lg[ck]:
                        cg = lg[ck][c_key]
                        cs = int(cg.attrs.get('start', 0))
                        ce = int(cg.attrs.get('end', 0))
                        cw = int(cg.attrs.get('width', 0))
                        if zt == 3:
                            ax.axvspan(cs, ce, alpha=0.15, color='red')
                            ax.annotate(f'{cw}bp', (cs, ax.get_ylim()[1] * 0.9),
                                        fontsize=6, color='red')

            ax.set_ylabel(f'Z-score\n({label})', fontsize=8)
            ax.set_xlim(0, ref_len)
            ax.axvspan(ps, pe, alpha=0.03, color='blue')

        axes[-1].set_xlabel('Position (bp)')
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig); plt.close(fig)

        # Separate page: Occupancy comparison (variant vs WT null mean)
        fig, axes = plt.subplots(len(bins), 1, figsize=(14, 3 * len(bins)), sharex=True)
        if len(bins) == 1:
            axes = [axes]
        fig.suptitle(f"Occupancy: Variant vs WT Null Mean — {vid_key}", fontsize=14, fontweight='bold')

        for i, label in enumerate(bins):
            s = sn(label)
            if s not in vg:
                continue
            lg = vg[s]
            ax = axes[i]
            vo = lg['variant_occupancy'][:]
            nm = lg['null_mean'][:]
            ns = lg['null_std'][:]
            positions = np.arange(len(vo))

            ax.fill_between(positions, nm - 2*ns, nm + 2*ns, alpha=0.15, color='steelblue',
                            label='WT null ±2σ')
            ax.plot(positions, nm, color='steelblue', lw=1, alpha=0.8, label='WT null mean')
            ax.plot(positions, vo, color='coral', lw=1.5, alpha=0.9, label='Variant')

            if vpos is not None:
                ax.axvline(vpos, color='red', ls='--', lw=1.5, alpha=0.6)

            ax.set_ylabel(f'Occupancy\n({label})', fontsize=8)
            ax.legend(fontsize=7, loc='upper right')
            ax.set_xlim(0, ref_len)

        axes[-1].set_xlabel('Position (bp)')
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Phase 3: Site-Level Test Results
# ============================================================================

def plot_site_tests(f, meta, pdf):
    if 'site_tests' not in f:
        return
    st = f['site_tests']
    bins = meta['bin_labels']

    for vid_key in st:
        vg = st[vid_key]
        # Collect all site results across bins
        all_sites = []
        for label in bins:
            s = sn(label)
            if s not in vg:
                continue
            for site_key in sorted(vg[s].keys()):
                sg = vg[s][site_key]
                attrs = dict(sg.attrs)
                if 'error' in attrs:
                    continue
                attrs['bin'] = label
                all_sites.append(attrs)

        if not all_sites:
            continue

        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        fig.suptitle(f"Phase 3: Site-Level Tests — {vid_key}", fontsize=14, fontweight='bold')

        # Forest plot: delta occupancy with OR
        ax = axes[0]
        labels_for_plot = []
        deltas = []
        pvals = []
        colors_plot = []
        for s in all_sites:
            lab = f"{s['bin']}\n{int(s['site_start'])}-{int(s['site_end'])}"
            labels_for_plot.append(lab)
            deltas.append(float(s.get('delta_occupancy', 0)))
            pvals.append(float(s.get('fisher_p', 1)))
            colors_plot.append(get_color(s['bin']))

        y = np.arange(len(labels_for_plot))
        ax.barh(y, deltas, color=colors_plot, alpha=0.7)
        ax.axvline(0, color='gray', ls='-', lw=0.5)
        ax.set_yticks(y); ax.set_yticklabels(labels_for_plot, fontsize=7)
        ax.set_xlabel('Δ Occupancy (Variant − WT)')
        ax.set_title('Effect Size')
        for i, p in enumerate(pvals):
            sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
            ax.annotate(f'p={p:.1e} {sig}', (deltas[i], i),
                        fontsize=6, va='center', ha='left' if deltas[i] >= 0 else 'right')

        # WT vs Variant occupancy
        ax = axes[1]
        wt_occs = [float(s.get('wt_occupancy', 0)) for s in all_sites]
        var_occs = [float(s.get('var_occupancy', 0)) for s in all_sites]
        ax.scatter(wt_occs, var_occs, c=colors_plot, s=80, alpha=0.8, edgecolors='black', lw=0.5)
        lim = max(max(wt_occs + [0.01]), max(var_occs + [0.01])) * 1.1
        ax.plot([0, lim], [0, lim], 'k--', alpha=0.3)
        ax.set_xlabel('WT Occupancy'); ax.set_ylabel('Variant Occupancy')
        ax.set_title('Occupancy: WT vs Variant')
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_aspect('equal')

        # P-value volcano
        ax = axes[2]
        for s in all_sites:
            d = float(s.get('delta_occupancy', 0))
            p = float(s.get('fisher_p', 1))
            log_p = -np.log10(max(p, 1e-300))
            c = get_color(s['bin'])
            ax.scatter(d, log_p, c=c, s=80, alpha=0.8, edgecolors='black', lw=0.5)
        ax.axhline(-np.log10(0.05), ls=':', color='gray', alpha=0.5, label='p=0.05')
        ax.axhline(-np.log10(0.01), ls=':', color='orange', alpha=0.5, label='p=0.01')
        ax.axvline(0, color='gray', ls='-', lw=0.5)
        ax.set_xlabel('Δ Occupancy'); ax.set_ylabel('-log10(p)')
        ax.set_title('Volcano Plot')
        ax.legend(fontsize=7)

        fig.tight_layout(rect=[0, 0, 1, 0.93])
        pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Phase 4: Distance Decay
# ============================================================================

def plot_distance_decay(f, meta, pdf):
    if 'distance_decay' not in f:
        return
    dd = f['distance_decay']
    bins = meta['bin_labels']

    for vid_key in dd:
        vg = dd[vid_key]
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle(f"Phase 4: Distance Decay — {vid_key}", fontsize=14, fontweight='bold')

        for label in bins:
            s = sn(label)
            if s not in vg:
                continue
            lg = vg[s]
            if 'distances' not in lg:
                continue
            dists = lg['distances'][:]
            mean_z = lg['mean_abs_z'][:]
            max_z = lg['max_abs_z'][:]
            n_pos = lg['n_positions'][:] if 'n_positions' in lg else np.ones_like(dists)
            mz = lg['mean_z'][:] if 'mean_z' in lg else np.zeros_like(dists)

            c = get_color(label)
            axes[0].plot(dists, mean_z, '-o', color=c, ms=3, label=label, alpha=0.8)
            axes[1].plot(dists, mz, '-o', color=c, ms=3, label=label, alpha=0.8)

        axes[0].set_xlabel('Distance from Variant (bp)')
        axes[0].set_ylabel('Mean |Z|')
        axes[0].set_title('Effect Magnitude vs Distance')
        axes[0].axhline(2, ls=':', color='gray', alpha=0.4)
        axes[0].legend(fontsize=8)

        axes[1].set_xlabel('Distance from Variant (bp)')
        axes[1].set_ylabel('Mean Z (signed)')
        axes[1].set_title('Effect Direction vs Distance')
        axes[1].axhline(0, ls='-', color='gray', alpha=0.3)
        axes[1].legend(fontsize=8)

        fig.tight_layout(rect=[0, 0, 1, 0.93])
        pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Summary Page
# ============================================================================

def plot_summary_page(f, meta, pdf):
    """Title page with key parameters and results."""
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.axis('off')

    lines = [
        f"MPRA Fiber-seq Calibration Report",
        f"",
        f"Version: {meta['version']}",
        f"Reference: {meta['ref_name']} ({meta['ref_length']} bp)",
        f"Promoter: {meta['prom_start']}-{meta['prom_end']} "
        f"({meta['prom_end'] - meta['prom_start']} bp)",
        f"Variant: {meta['variant_id']} at position {meta['variant_position']}",
        f"Total reads parsed: {meta['parse_parsed']:,} / {meta['parse_total']:,}",
        f"Bins: {', '.join(meta['bin_labels'])}",
        f"",
    ]

    # Add calibration recommendations if available
    if 'nucleosome_calibration' in f:
        nc = f['nucleosome_calibration']
        lines.append(f"Nucleosome count KS test: stat={nc.attrs.get('ks_stat',0):.4f}, "
                     f"p={nc.attrs.get('ks_pval',1):.2e}")

    if 'wt_null_analysis' in f:
        wn = f['wt_null_analysis']
        for label in meta['bin_labels']:
            s = sn(label)
            if s in wn:
                fz2 = float(wn[s].attrs.get('frac_z_gt_2', 0))
                c95 = float(wn[s].attrs.get('cluster_size_95th', 0))
                lines.append(f"  {label}: null |Z|>2 rate = {fz2:.3%}, "
                             f"cluster 95th = {c95:.0f} bp")

    if 'binding_sites' in f:
        bs = f['binding_sites']
        lines.append(f"")
        lines.append(f"Binding sites called:")
        for label in meta['bin_labels']:
            s = sn(label)
            if s in bs:
                n_sites = len(bs[s])
                lines.append(f"  {label}: {n_sites} sites")

    text = '\n'.join(lines)
    ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    fig.tight_layout()
    pdf.savefig(fig); plt.close(fig)


# ============================================================================
# Main
# ============================================================================

def main():
    p = argparse.ArgumentParser(description='Visualize MPRA Fiber-seq calibration results')
    p.add_argument('h5_file', help='HDF5 output from mpra_fiberseq_analysis.py')
    p.add_argument('--output', '-o', default=None,
                   help='Output PDF (default: <input>_plots.pdf)')
    args = p.parse_args()

    h5_path = Path(args.h5_file)
    if args.output:
        pdf_path = Path(args.output)
    else:
        pdf_path = h5_path.with_name(h5_path.stem + '_plots.pdf')

    print(f"Reading: {h5_path}")
    print(f"Output:  {pdf_path}")

    with h5py.File(h5_path, 'r') as f:
        meta = load_meta(f)
        print(f"Reference: {meta['ref_name']} ({meta['ref_length']} bp)")
        print(f"Bins: {meta['bin_labels']}")
        print(f"Variant: {meta['variant_id']} at {meta['variant_position']}")
        print(f"Groups in HDF5: {list(f.keys())}")

        with PdfPages(str(pdf_path)) as pdf:
            plot_summary_page(f, meta, pdf)
            plot_nuc_calibration(f, meta, pdf)
            plot_quality_calibration(f, meta, pdf)
            plot_power_analysis(f, meta, pdf)
            plot_wt_null(f, meta, pdf)
            plot_sequence_context(f, meta, pdf)
            plot_wt_landscape(f, meta, pdf)
            plot_position_tests(f, meta, pdf)
            plot_site_tests(f, meta, pdf)
            plot_distance_decay(f, meta, pdf)

    print(f"\nDone! {pdf_path}")


if __name__ == '__main__':
    main()
