#!/usr/bin/env python3
"""
mpra_browser.py — Interactive MPRA Fiber-seq Library Browser

Dash/Plotly web app for browsing saturation mutagenesis results.
Launch: python mpra_browser.py --h5 library.h5 [--mpra mpra_data.tsv] [--port 8050]

Tabs:
  1. WT Landscape    — Baseline WT protein footprint occupancy
  2. Variant Browser — Step through variants with occupancy, Δ, significance
  3. Library Heatmap — All variants × positions heatmap
  4. MPRA Overlay    — Functional MPRA scores + footprint integration
  5. Comparison      — Side-by-side variant comparison

Author: Ben / Stergachis Lab, University of Washington
"""

import argparse
import re
import sys
import warnings

import h5py
import numpy as np
import pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

from dash import Dash, html, dcc, callback, Input, Output, State, no_update

warnings.filterwarnings('ignore')

def _safe_decode(val):
    """Decode bytes to str, pass through strings unchanged."""
    if isinstance(val, bytes):
        return val.decode()
    return str(val)


def _hex_to_rgba(hex_color, alpha=0.2):
    """Convert hex color to rgba string for Plotly."""
    hex_color = hex_color.lstrip('#')
    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    return f'rgba({r},{g},{b},{alpha})'


# ============================================================================
# Argument Parsing
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='MPRA Fiber-seq Library Browser')
    p.add_argument('--h5', required=True, help='Library HDF5 file')
    p.add_argument('--mpra', default=None, help='MPRA functional data TSV/CSV')
    p.add_argument('--port', type=int, default=8050)
    p.add_argument('--host', default='127.0.0.1')
    p.add_argument('--debug', action='store_true')

    # Coordinate mapping
    p.add_argument('--chrom', default=None, help='Genomic chromosome (e.g. chr19)')
    p.add_argument('--genomic-start', type=int, default=None,
                   help='Genomic start (1-based inclusive)')
    p.add_argument('--genomic-end', type=int, default=None,
                   help='Genomic end (1-based inclusive)')
    p.add_argument('--strand', default='-', choices=['+', '-'],
                   help='Strand orientation (default: -)')
    return p.parse_args()


# ============================================================================
# Data Loading
# ============================================================================

class LibraryData:
    """Loads and manages all data from the HDF5 + optional MPRA file."""

    BIN_COLORS = {
        'sub-TF_10-19bp': '#66c2a5', 'TF_20-40bp': '#fc8d62',
        'PIC_41-80bp': '#8da0cb', 'Nucleosome_81plusbp': '#e78ac3',
    }
    BIN_SHORT = {
        'sub-TF_10-19bp': 'sub-TF', 'TF_20-40bp': 'TF',
        'PIC_41-80bp': 'PIC', 'Nucleosome_81plusbp': 'NUC',
    }

    def __init__(self, h5_path, mpra_path=None,
                 chrom=None, genomic_start=None, genomic_end=None,
                 strand='-'):
        self.h5 = h5py.File(h5_path, 'r')
        self._load_metadata()
        self._load_wt_occupancy()
        self._load_summary_table()
        self._setup_coordinate_mapping(chrom, genomic_start, genomic_end,
                                        strand)
        self._build_variant_display_info()
        if mpra_path:
            self._load_mpra(mpra_path)
            self._diagnose_mpra_matching()
        else:
            self.mpra_df = None

    def _load_metadata(self):
        meta = self.h5['metadata']
        raw_labels = meta.attrs['bin_labels']
        self.bin_labels = [b.decode() if isinstance(b, bytes) else str(b)
                           for b in raw_labels]
        self.prom_s = int(meta.attrs['promoter_start'])
        self.prom_e = int(meta.attrs['promoter_end'])
        self.a_start = int(meta.attrs.get('analysis_start', self.prom_s))
        self.a_end = int(meta.attrs.get('analysis_end', self.prom_e))
        self.analysis_length = self.a_end - self.a_start
        self.ref_name = meta.attrs.get('reference_name', '')
        self.ref_length = int(meta.attrs.get('reference_length', 0))
        self.n_variants = int(meta.attrs.get('n_variants_tested', 0))
        print(f"Loaded: {self.n_variants} variants, "
              f"{self.analysis_length} bp analysis region")

    def _load_wt_occupancy(self):
        self.wt_occ = {}
        # Try top-level wt_occupancy first (new format)
        if 'wt_occupancy' in self.h5:
            wt_grp = self.h5['wt_occupancy']
            for label in self.bin_labels:
                # Try raw name, then sanitized
                for key in [label, self._safe_name(label)]:
                    if key in wt_grp:
                        self.wt_occ[label] = wt_grp[key][:]
                        break
        else:
            # Fallback: grab from first null calibration depth
            nc = self.h5['null_calibration']
            first_depth = list(nc.keys())[0]
            depth_grp = nc[first_depth]
            for label in self.bin_labels:
                for key in [label, self._safe_name(label)]:
                    if key in depth_grp:
                        subgrp = depth_grp[key]
                        if 'wt_occ' in subgrp:
                            self.wt_occ[label] = subgrp['wt_occ'][:]
                            break

    def _load_summary_table(self):
        sg = self.h5['summary']
        self.summary = pd.DataFrame({
            'variant_id': [_safe_decode(v) for v in sg['variant_ids'][:]],
            'position': sg['positions'][:],
            'n_reads': sg['n_reads'][:],
            'best_cluster_p': sg['best_cluster_p'][:],
            'fdr_q': sg['variant_fdr_q'][:],
        })
        # Add ref/alt if available
        if 'ref_bases' in sg:
            self.summary['ref'] = [_safe_decode(v) for v in sg['ref_bases'][:]]
            self.summary['alt'] = [_safe_decode(v) for v in sg['alt_bases'][:]]
            self.summary['change_type'] = [_safe_decode(v) for v in sg['change_types'][:]]
        if 'n_reads_raw' in sg:
            self.summary['n_reads_raw'] = sg['n_reads_raw'][:]

        # ── Build hdf5_key mapping by scanning actual /variants/ groups ──
        vg = self.h5['variants']
        actual_keys = list(vg.keys())

        # Map variant_id -> actual HDF5 group key
        vid_to_key = {}
        for gkey in actual_keys:
            vid = _safe_decode(vg[gkey].attrs.get('variant_id', ''))
            if vid:
                vid_to_key[vid] = gkey
        self.summary['hdf5_key'] = self.summary['variant_id'].map(
            vid_to_key).fillna('')

        # ── Discover bin label -> actual subgroup name mapping ──
        self._bin_key_map = {}
        if actual_keys:
            first_v = vg[actual_keys[0]]
            subgroups = list(first_v.keys())
            for label in self.bin_labels:
                if label in subgroups:
                    self._bin_key_map[label] = label
                else:
                    sl = self._safe_name(label)
                    if sl in subgroups:
                        self._bin_key_map[label] = sl
                    else:
                        # Fuzzy match
                        for sg_name in subgroups:
                            if (label.split('_')[0] in sg_name or
                                    sg_name.split('_')[0] in label):
                                self._bin_key_map[label] = sg_name
                                break
            print(f"  Bin key map: {self._bin_key_map}")
        else:
            self._bin_key_map = {l: l for l in self.bin_labels}

        # Per-bin summary stats
        for label in self.bin_labels:
            sl = self._safe_name(label)
            short = self.BIN_SHORT.get(label, label)
            if f'{sl}_max_abs_delta' in sg:
                self.summary[f'{short}_max_delta'] = sg[f'{sl}_max_abs_delta'][:]
            if f'{sl}_n_sig_clusters' in sg:
                self.summary[f'{short}_n_clusters'] = sg[f'{sl}_n_sig_clusters'][:]

        # Load cluster table if available
        if 'clusters' in self.h5:
            cg = self.h5['clusters']
            self.clusters = pd.DataFrame({
                'variant_idx': cg['variant_idx'][:],
                'variant_id': [_safe_decode(v) for v in cg['variant_ids'][:]],
                'bin_label': [_safe_decode(v) for v in cg['bin_labels'][:]],
                'abs_start': cg['abs_start'][:],
                'abs_end': cg['abs_end'][:],
                'width': cg['width'][:],
                'sum_abs_delta': cg['sum_abs_delta'][:],
                'mean_signed_delta': cg['mean_signed_delta'][:],
                'direction': [_safe_decode(v) for v in cg['direction'][:]],
                'sum_p': cg['sum_p'][:],
            })
        else:
            self.clusters = pd.DataFrame()

        n_matched = (self.summary['hdf5_key'] != '').sum()
        print(f"  Summary: {len(self.summary)} variants, "
              f"{n_matched} matched to HDF5 groups")
        if len(self.clusters) > 0:
            print(f"  Clusters: {len(self.clusters)} total")

    def _setup_coordinate_mapping(self, chrom, genomic_start, genomic_end,
                                   strand):
        self.chrom = chrom
        self.genomic_start = genomic_start
        self.genomic_end = genomic_end
        self.strand = strand
        self.use_genomic = (chrom is not None and genomic_start is not None
                            and genomic_end is not None
                            and genomic_end > genomic_start)

        if self.use_genomic:
            prom_span = self.prom_e - self.prom_s
            gen_span = genomic_end - genomic_start + 1
            if prom_span != gen_span:
                print(f"WARNING: promoter span ({prom_span}) != genomic span "
                      f"({gen_span}). Disabling coordinate mapping.")
                self.use_genomic = False
            else:
                print(f"  Coordinate mapping: {chrom}:{genomic_start:,}-"
                      f"{genomic_end:,} ({strand} strand)")

    def _build_variant_display_info(self):
        """Add display labels and positions to summary table."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                      'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}

        display_labels = []
        display_positions = []
        for _, row in self.summary.iterrows():
            vid = row['variant_id']
            plasmid_pos = row['position']
            if self.use_genomic and plasmid_pos >= 0:
                dpos = self._plasmid_to_display(plasmid_pos)
                # Convert label
                match = re.match(r'(\d+):(.*)', vid)
                if match and self.strand == '-':
                    change = match.group(2)
                    conv = ''.join(complement.get(c, c) for c in change)
                    display_labels.append(f'{dpos}:{conv}')
                elif match:
                    display_labels.append(f'{dpos}:{match.group(2)}')
                else:
                    display_labels.append(vid)
                display_positions.append(dpos)
            else:
                display_labels.append(vid)
                display_positions.append(plasmid_pos)

        self.summary['display_label'] = display_labels
        self.summary['display_pos'] = display_positions

    def _plasmid_to_display(self, pos_0based):
        if not self.use_genomic:
            return pos_0based
        offset = pos_0based - self.prom_s
        if self.strand == '+':
            return self.genomic_start + offset
        else:
            return self.genomic_end - offset

    def get_display_positions(self):
        if not self.use_genomic:
            return np.arange(self.a_start, self.a_end)
        return np.arange(self.genomic_start, self.genomic_end + 1)

    def reorient(self, arr):
        if not self.use_genomic or self.strand == '+':
            return arr
        return arr[::-1]

    def get_xlabel(self):
        if not self.use_genomic:
            return 'Plasmid Position (bp)'
        return f'{self.chrom} ({self.strand} strand, hg38)'

    def get_variant_data(self, hdf5_key):
        """Load per-variant arrays for one variant."""
        vg = self.h5['variants']
        if not hdf5_key or hdf5_key not in vg:
            return None
        v = vg[hdf5_key]
        data = {
            'variant_id': _safe_decode(v.attrs.get('variant_id', '')),
            'n_reads': int(v.attrs.get('n_nc_matched', 0)),
            'best_p': float(v.attrs.get('best_cluster_p', 1.0)),
            'fdr_q': float(v.attrs.get('variant_fdr_q', 1.0)),
        }
        for label in self.bin_labels:
            # Use discovered bin key mapping
            bin_key = self._bin_key_map.get(label, label)
            if bin_key not in v:
                continue
            lg = v[bin_key]
            data[label] = {
                'delta_obs': self.reorient(lg['delta_obs'][:]),
                'variant_occ': self.reorient(lg['variant_occ'][:]),
                'empirical_p': self.reorient(lg['empirical_p'][:]),
                'q_values': self.reorient(lg['q_values'][:]),
                'z_scores': self.reorient(lg['z_scores'][:]),
                'n_sig_positions': int(lg.attrs.get(
                    'n_sig_positions_fdr10', 0)),
                'max_abs_delta': float(lg.attrs.get('max_abs_delta', 0)),
            }
            # Load clusters
            sc = lg.get('significant_clusters')
            clusters = []
            if sc is not None:
                for cname in sorted(sc.keys()):
                    c = sc[cname]
                    cs = self._plasmid_to_display(
                        int(c.attrs['abs_start']))
                    ce = self._plasmid_to_display(
                        int(c.attrs['abs_end']))
                    clusters.append({
                        'start': min(cs, ce),
                        'end': max(cs, ce),
                        'width': int(c.attrs['width']),
                        'sum_abs_delta': float(
                            c.attrs['sum_abs_delta']),
                        'sum_p': float(c.attrs.get('sum_p', 1.0)),
                        'direction': _safe_decode(
                            c.attrs.get('direction', '')),
                    })
            data[label]['clusters'] = clusters
        return data

    def _load_mpra(self, mpra_path):
        """Load MPRA functional data from TSV/CSV."""
        sep = ',' if mpra_path.endswith('.csv') else '\t'
        try:
            self.mpra_df = pd.read_csv(mpra_path, sep=sep)
            # Standardize column names
            col_map = {}
            for col in self.mpra_df.columns:
                cl = col.lower().strip()
                if cl == 'position':
                    col_map[col] = 'position'
                elif cl == 'chromosome' or cl == 'chrom':
                    col_map[col] = 'chrom'
                elif cl == 'ref':
                    col_map[col] = 'ref'
                elif cl == 'alt':
                    col_map[col] = 'alt'
                elif cl == 'value':
                    col_map[col] = 'value'
                elif cl == 'p-value' or cl == 'pvalue' or cl == 'p_value':
                    col_map[col] = 'pvalue'
                elif cl == 'tags':
                    col_map[col] = 'tags'
                elif cl == 'dna':
                    col_map[col] = 'dna_count'
                elif cl == 'rna':
                    col_map[col] = 'rna_count'
            self.mpra_df = self.mpra_df.rename(columns=col_map)

            # Compute significance
            if 'pvalue' in self.mpra_df.columns:
                self.mpra_df['significant'] = self.mpra_df['pvalue'] < 0.05

            print(f"  MPRA data: {len(self.mpra_df)} variants loaded")
            print(f"  Columns: {list(self.mpra_df.columns)}")
        except Exception as e:
            print(f"  WARNING: Could not load MPRA data: {e}")
            self.mpra_df = None

    def _diagnose_mpra_matching(self):
        """Print diagnostic info about MPRA-to-footprint matching."""
        if self.mpra_df is None:
            return
        n_matched = 0
        n_exact = 0
        n_pos_only = 0
        n_unmatched = 0
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                      'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
        for _, row in self.summary.iterrows():
            m = _get_mpra_for_variant(self, row)
            if m is not None:
                n_matched += 1
                # Check if it was an exact ref/alt match
                ref = row.get('ref', '')
                alt = row.get('alt', '')
                if self.use_genomic and self.strand == '-' and ref and alt:
                    ref = complement.get(ref, ref)
                    alt = complement.get(alt, alt)
                if (ref and alt and 'ref' in self.mpra_df.columns and
                        m.get('ref', '').upper() == ref.upper() and
                        m.get('alt', '').upper() == alt.upper()):
                    n_exact += 1
                else:
                    n_pos_only += 1
            else:
                n_unmatched += 1

        total = len(self.summary)
        print(f"  MPRA matching: {n_matched}/{total} variants matched "
              f"({n_exact} exact ref/alt, {n_pos_only} position-only, "
              f"{n_unmatched} unmatched)")
        if n_unmatched > total * 0.5:
            print(f"  WARNING: >50% unmatched — check coordinate system!")
            # Show some examples
            for _, row in self.summary.head(5).iterrows():
                m = _get_mpra_for_variant(self, row)
                status = 'MATCHED' if m is not None else 'NO MATCH'
                print(f"    {row['display_label']} (pos={row['display_pos']}) "
                      f"-> {status}")

    @staticmethod
    def _safe_name(name):
        return (name.replace(' ', '_').replace('(', '').replace(')', '')
                    .replace('+', 'plus').replace('>', 'to')
                    .replace(':', '_'))

    @staticmethod
    def _safe_variant_name(vid):
        return (vid.replace(':', '_').replace('>', 'to')
                    .replace('+', 'plus').replace(' ', '')
                    .replace('"', '').replace('[', '')
                    .replace(']', '').replace(',', '_')
                    .replace('/', '_'))



# ============================================================================
# Helpers
# ============================================================================

# Colorscale options for heatmap
COLORSCALE_OPTIONS = [
    {'label': 'RdBu (diverging)', 'value': 'RdBu_r'},
    {'label': 'Hot', 'value': 'Hot_r'},
    {'label': 'Viridis', 'value': 'Viridis'},
    {'label': 'Plasma', 'value': 'Plasma'},
    {'label': 'YlOrRd', 'value': 'YlOrRd'},
    {'label': 'PiYG (diverging)', 'value': 'PiYG'},
    {'label': 'BrBG (diverging)', 'value': 'BrBG'},
    {'label': 'Greys', 'value': 'Greys'},
]

def _get_mpra_for_variant(data, variant_row):
    """Find MPRA data matching a footprint variant.

    Matches on display_pos (genomic if mapping active, else plasmid).
    For minus strand, complements the footprint ref/alt before comparing
    to MPRA ref/alt (which are in genomic orientation).
    """
    if data.mpra_df is None or 'position' not in data.mpra_df.columns:
        return None
    display_pos = variant_row.get('display_pos', 0)
    if display_pos == 0:
        return None
    matches = data.mpra_df[data.mpra_df['position'] == display_pos]
    if len(matches) == 0:
        return None

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}

    ref = variant_row.get('ref', '')
    alt = variant_row.get('alt', '')

    # On minus strand, complement the plasmid bases to match genomic MPRA bases
    if data.use_genomic and data.strand == '-' and ref and alt:
        ref = complement.get(ref, ref)
        alt = complement.get(alt, alt)

    if ref and alt and 'ref' in matches.columns and 'alt' in matches.columns:
        exact = matches[(matches['ref'].str.upper() == ref.upper()) &
                        (matches['alt'].str.upper() == alt.upper())]
        if len(exact) > 0:
            return exact.iloc[0]
    # Fallback: return first match at position
    return matches.iloc[0]


# ============================================================================
# Tab 1: WT Landscape (MPRA on top, crosshair on hover)
# ============================================================================

def build_wt_landscape(data, show_mpra=True, hover_pos=None):
    """WT landscape with MPRA on top row."""
    x = data.get_display_positions()
    n_bins = len(data.bin_labels)
    has_mpra = show_mpra and data.mpra_df is not None and 'value' in data.mpra_df.columns

    n_rows = (1 if has_mpra else 0) + n_bins
    row_heights = ([0.6] if has_mpra else []) + [1.0] * n_bins
    titles = (['MPRA Log₂ Effect'] if has_mpra else []) + \
             [data.BIN_SHORT.get(l, l) for l in data.bin_labels]

    fig = make_subplots(
        rows=n_rows, cols=1, shared_xaxes=True,
        vertical_spacing=0.03, row_heights=row_heights,
        subplot_titles=titles)

    row_offset = 1 if has_mpra else 0

    # MPRA on top
    if has_mpra:
        df = data.mpra_df
        is_sig = df.get('significant', pd.Series([True] * len(df)))
        colors_mpra = np.where(is_sig, '#2ca02c', '#cccccc')
        fig.add_trace(go.Bar(
            x=df['position'], y=df['value'],
            marker_color=colors_mpra.tolist(),
            marker_line_width=0, name='MPRA',
            hovertemplate='pos: %{x:,}<br>log₂: %{y:.2f}<extra></extra>',
        ), row=1, col=1)
        fig.add_hline(y=0, line=dict(color='black', width=0.5), row=1, col=1)
        fig.update_yaxes(title_text='Log₂ effect', row=1, col=1)

    # Footprint bins
    for i, label in enumerate(data.bin_labels):
        r = i + 1 + row_offset
        wt = data.reorient(data.wt_occ.get(label, np.zeros(data.analysis_length)))
        color = data.BIN_COLORS.get(label, 'gray')
        fig.add_trace(go.Scatter(
            x=x, y=wt, mode='lines',
            fill='tozeroy', fillcolor=_hex_to_rgba(color, 0.2),
            line=dict(color=color, width=1.2),
            name=data.BIN_SHORT.get(label, label),
            hovertemplate='pos: %{x:,}<br>occ: %{y:.4f}<extra></extra>',
        ), row=r, col=1)
        fig.update_yaxes(title_text='Occupancy', row=r, col=1)

    # Hover crosshair line across all subplots
    if hover_pos is not None:
        for r in range(1, n_rows + 1):
            fig.add_vline(x=hover_pos,
                          line=dict(color='red', width=1, dash='dot'),
                          opacity=0.6, row=r, col=1)

    fig.update_xaxes(title_text=data.get_xlabel(), row=n_rows, col=1)
    fig.update_layout(
        height=200 * n_rows,
        title_text='WT Protein Footprint Landscape',
        showlegend=False,
        margin=dict(l=60, r=20, t=60, b=40),
        hovermode='x unified',
    )
    return fig


# ============================================================================
# Tab 2: Variant Browser
# ============================================================================

def build_variant_plot(data, variant_idx, delta_mode='gray'):
    """Per-variant plots with variant position dotted line on ALL plots."""
    row_info = data.summary.iloc[variant_idx]
    hkey = row_info['hdf5_key']
    vdata = data.get_variant_data(hkey)

    if vdata is None:
        fig = go.Figure()
        fig.add_annotation(text=f"No data for variant: {row_info['display_label']}",
                           xref="paper", yref="paper", x=0.5, y=0.5,
                           showarrow=False, font=dict(size=16))
        return fig, ""

    x = data.get_display_positions()
    n_bins = len(data.bin_labels)
    display_label = row_info['display_label']
    vpos = row_info['display_pos']

    mpra_match = _get_mpra_for_variant(data, row_info)

    fig = make_subplots(
        rows=n_bins, cols=3, shared_xaxes=True,
        horizontal_spacing=0.07, vertical_spacing=0.06,
        column_widths=[0.38, 0.34, 0.28])

    cluster_text_parts = []
    all_deltas_max = []

    for i, label in enumerate(data.bin_labels):
        color = data.BIN_COLORS.get(label, 'gray')
        short = data.BIN_SHORT.get(label, label)
        r = i + 1

        if label not in vdata:
            continue

        ld = vdata[label]
        wt_occ = data.reorient(data.wt_occ.get(label,
                               np.zeros(data.analysis_length)))
        delta = ld['delta_obs']
        var_occ = ld['variant_occ']
        emp_p = ld['empirical_p']
        q_vals = ld['q_values']
        clusters = ld.get('clusters', [])
        all_deltas_max.append(np.max(np.abs(delta)))

        # Col 1: Occupancy
        fig.add_trace(go.Scatter(
            x=x, y=wt_occ, mode='lines',
            line=dict(color='#4a86c8', width=1.2),
            name='WT', legendgroup='wt', showlegend=(i == 0),
            hovertemplate='%{x:,}: %{y:.4f}<extra>WT</extra>',
        ), row=r, col=1)
        fig.add_trace(go.Scatter(
            x=x, y=var_occ, mode='lines',
            line=dict(color='#e85d50', width=1.2),
            name='Variant', legendgroup='var', showlegend=(i == 0),
            hovertemplate='%{x:,}: %{y:.4f}<extra>Variant</extra>',
        ), row=r, col=1)
        for cl in clusters:
            fig.add_vrect(x0=cl['start'], x1=cl['end'],
                          fillcolor='red', opacity=0.1,
                          layer='below', line_width=0, row=r, col=1)
        fig.update_yaxes(title_text=f'{short} Occ', row=r, col=1,
                         title_font_size=10)

        # Col 2: Δ
        sig_mask = q_vals < 0.10
        if delta_mode == 'gray':
            fig.add_trace(go.Scatter(
                x=x, y=delta, mode='lines',
                line=dict(color='lightgray', width=0.8),
                showlegend=False,
                hovertemplate='%{x:,}: Δ=%{y:.4f}<extra>NS</extra>',
            ), row=r, col=2)
            delta_sig = np.where(sig_mask, delta, np.nan)
            fig.add_trace(go.Scatter(
                x=x, y=delta_sig, mode='lines',
                line=dict(color=color, width=1.5),
                fill='tozeroy', fillcolor=_hex_to_rgba(color, 0.2),
                showlegend=False,
                hovertemplate='%{x:,}: Δ=%{y:.4f}<extra>FDR<0.10</extra>',
            ), row=r, col=2)
        elif delta_mode == 'blank':
            delta_blank = np.where(sig_mask, delta, 0.0)
            fig.add_trace(go.Scatter(
                x=x, y=delta_blank, mode='lines',
                line=dict(color=color, width=1.2),
                fill='tozeroy', fillcolor=_hex_to_rgba(color, 0.2),
                showlegend=False,
            ), row=r, col=2)
        else:
            fig.add_trace(go.Scatter(
                x=x, y=delta, mode='lines',
                line=dict(color=color, width=1),
                fill='tozeroy', fillcolor=_hex_to_rgba(color, 0.2),
                showlegend=False,
            ), row=r, col=2)
        for cl in clusters:
            fig.add_vrect(x0=cl['start'], x1=cl['end'],
                          fillcolor='red', opacity=0.1,
                          layer='below', line_width=0, row=r, col=2)
        fig.add_hline(y=0, line=dict(color='gray', width=0.5), row=r, col=2)
        fig.update_yaxes(title_text=f'{short} Δ', row=r, col=2,
                         title_font_size=10, title_standoff=5)

        # Col 3: -log10(p)
        neg_log_p = -np.log10(np.maximum(emp_p, 1e-10))
        fig.add_trace(go.Scatter(
            x=x, y=neg_log_p, mode='lines',
            line=dict(color=color, width=0.8), showlegend=False,
            hovertemplate='%{x:,}: -log₁₀(p)=%{y:.2f}<extra></extra>',
        ), row=r, col=3)
        fig.add_hline(y=-np.log10(0.05),
                      line=dict(color='orange', width=0.5, dash='dash'),
                      row=r, col=3)
        fig.add_hline(y=-np.log10(0.01),
                      line=dict(color='red', width=0.5, dash='dash'),
                      row=r, col=3)
        fig.update_yaxes(title_text='-log₁₀(p)', row=r, col=3,
                         title_font_size=10, title_standoff=5)

        for cl in clusters:
            cluster_text_parts.append(
                f"**{short}**: [{cl['start']:,}–{cl['end']:,}] "
                f"w={cl['width']}bp Σ|Δ|={cl['sum_abs_delta']:.3f} "
                f"p={cl['sum_p']:.4f} dir={cl['direction']}")

    # Variant position dotted line on ALL subplots using explicit shapes
    if vpos is not None and vpos != 0:
        # Plotly subplot axes: row 1 col 1 = x1/y1, row 1 col 2 = x2/y2, etc.
        # For n_bins rows × 3 cols: axis index = (row-1)*3 + col
        for r in range(1, n_bins + 1):
            for c in range(1, 4):
                ax_idx = (r - 1) * 3 + c
                xref = 'x' if ax_idx == 1 else f'x{ax_idx}'
                yref = 'y' if ax_idx == 1 else f'y{ax_idx}'
                fig.add_shape(
                    type='line', x0=vpos, x1=vpos, y0=0, y1=1,
                    xref=xref, yref=f'{yref} domain',
                    line=dict(color='red', width=1.5, dash='dot'),
                    opacity=0.6,
                )

    # Standardize Δ y-axis
    if all_deltas_max:
        max_delta = max(all_deltas_max) * 1.1
        for i in range(n_bins):
            fig.update_yaxes(range=[-max_delta, max_delta],
                             row=i + 1, col=2)

    for col in [1, 2, 3]:
        fig.update_xaxes(title_text=data.get_xlabel(), row=n_bins, col=col)

    sig_marker = ('★★★' if row_info['fdr_q'] < 0.05
                  else ('★★' if row_info['fdr_q'] < 0.10
                        else ('★' if row_info['fdr_q'] < 0.20 else '')))
    title = (f"{display_label}  |  n={row_info['n_reads']}  |  "
             f"p={row_info['best_cluster_p']:.4f}  |  "
             f"FDR q={row_info['fdr_q']:.4f} {sig_marker}")
    if mpra_match is not None:
        mpra_val = mpra_match.get('value', float('nan'))
        mpra_sig = '(sig)' if mpra_match.get('significant', False) else '(NS)'
        title += f"  |  MPRA log₂={mpra_val:.2f} {mpra_sig}"

    fig.update_layout(
        height=220 * n_bins, title_text=title, title_font_size=13,
        showlegend=True,
        legend=dict(orientation='h', yanchor='bottom', y=1.02,
                    xanchor='right', x=1),
        margin=dict(l=50, r=20, t=80, b=40),
    )
    cluster_md = "\n\n".join(cluster_text_parts) if cluster_text_parts else "*No significant clusters*"
    return fig, cluster_md


# ============================================================================
# Tab 3: Library Heatmap (with direction filter, scale bar, colorscale)
# ============================================================================

def build_heatmap(data, bin_label, metric='-log10(p)',
                  direction_filter='all', zmin_override=None,
                  zmax_override=None, colorscale_name=None,
                  highlight_variant_idx=None):
    """Library heatmap with filtering and interactive controls."""
    bin_key = data._bin_key_map.get(bin_label, bin_label)
    short = data.BIN_SHORT.get(bin_label, bin_label)
    x = data.get_display_positions()

    sorted_df = data.summary.sort_values('display_pos', ascending=True).reset_index(drop=True)
    n_var = len(sorted_df)
    matrix = np.zeros((n_var, data.analysis_length), dtype=np.float32)

    for row_i, (_, row) in enumerate(sorted_df.iterrows()):
        hkey = row['hdf5_key']
        vg = data.h5['variants']
        if not hkey or hkey not in vg:
            continue
        v = vg[hkey]
        if bin_key not in v:
            continue
        lg = v[bin_key]

        if metric == '-log10(p)':
            raw = lg['empirical_p'][:]
            vals = -np.log10(np.maximum(raw, 1e-10))
        elif metric == 'Δ (signed)':
            vals = lg['delta_obs'][:]
        elif metric == '|Δ|':
            vals = np.abs(lg['delta_obs'][:])
        elif metric == 'Z-score':
            vals = lg['z_scores'][:]
        else:
            vals = lg['delta_obs'][:]

        vals = data.reorient(vals)

        # Direction filter
        if direction_filter == 'decrease':
            vals = np.where(data.reorient(lg['delta_obs'][:]) < 0, vals, 0.0)
        elif direction_filter == 'increase':
            vals = np.where(data.reorient(lg['delta_obs'][:]) > 0, vals, 0.0)

        matrix[row_i, :] = vals

    # Y-axis labels with base change + MPRA
    y_labels = []
    has_mpra = data.mpra_df is not None and 'value' in data.mpra_df.columns
    for _, row in sorted_df.iterrows():
        ref = row.get('ref', '')
        alt = row.get('alt', '')
        ct = row.get('change_type', '')
        if ref and alt:
            lbl = f"{row['display_pos']}:{ref}>{alt}"
        elif ct == 'del':
            lbl = f"{row['display_pos']}:Δ{ref}"
        elif ct == 'ins':
            lbl = f"{row['display_pos']}:+{alt}"
        else:
            lbl = row['display_label']
        if has_mpra:
            m = _get_mpra_for_variant(data, row)
            if m is not None:
                lbl += f" [{m.get('value', 0):+.1f}]"
        y_labels.append(lbl)

    # Colorscale and range
    if colorscale_name is None:
        if metric in ('Δ (signed)', 'Z-score'):
            colorscale_name = 'RdBu_r'
        elif metric == '-log10(p)':
            colorscale_name = 'Hot_r'
        else:
            colorscale_name = 'YlOrRd'

    nonzero = matrix[matrix != 0]
    if len(nonzero) == 0:
        nonzero = np.array([0, 1])

    if metric in ('Δ (signed)', 'Z-score'):
        vm = np.percentile(np.abs(nonzero), 99)
        default_zmin, default_zmax = -vm, vm
    elif metric == '-log10(p)':
        default_zmin = 0
        default_zmax = min(np.percentile(nonzero, 99), 4)
    else:
        default_zmin = 0
        default_zmax = np.percentile(nonzero, 99)

    zmin = zmin_override if zmin_override is not None else default_zmin
    zmax = zmax_override if zmax_override is not None else default_zmax

    # Build figure: MPRA sidebar + heatmap
    if has_mpra:
        fig = make_subplots(
            rows=1, cols=2, shared_yaxes=True,
            column_widths=[0.08, 0.92],
            horizontal_spacing=0.01)

        mpra_vals = []
        mpra_hover = []
        for _, row in sorted_df.iterrows():
            m = _get_mpra_for_variant(data, row)
            if m is not None:
                val = m.get('value', 0.0)
                mpra_vals.append(val)
                mpra_ref = m.get('ref', '?')
                mpra_alt = m.get('alt', '?')
                mpra_pos = m.get('position', '?')
                sig_str = 'sig' if m.get('significant', False) else 'NS'
                mpra_hover.append(
                    f"Footprint: {row['display_label']}<br>"
                    f"MPRA: {mpra_pos}:{mpra_ref}>{mpra_alt}<br>"
                    f"log₂={val:.2f} ({sig_str})")
            else:
                mpra_vals.append(0.0)
                mpra_hover.append(
                    f"Footprint: {row['display_label']}<br>"
                    f"No MPRA match")

        mpra_colors = ['#d62728' if v < -0.5 else '#2ca02c' if v > 0.5
                       else '#cccccc' for v in mpra_vals]
        fig.add_trace(go.Bar(
            y=list(range(n_var)), x=mpra_vals, orientation='h',
            marker_color=mpra_colors, marker_line_width=0,
            customdata=mpra_hover,
            hovertemplate='%{customdata}<extra></extra>',
            showlegend=False,
        ), row=1, col=1)
        fig.update_xaxes(title_text='MPRA', row=1, col=1, autorange='reversed')

        fig.add_trace(go.Heatmap(
            z=matrix, x=x, y=list(range(n_var)),
            colorscale=colorscale_name, zmin=zmin, zmax=zmax,
            hovertemplate=('pos: %{x:,}<br>%{customdata}<br>'
                           f'{metric}: %{{z:.3f}}<extra></extra>'),
            customdata=np.array([[y_labels[r]] * data.analysis_length
                                 for r in range(n_var)]),
            colorbar=dict(title=metric),
        ), row=1, col=2)
        fig.update_xaxes(title_text=data.get_xlabel(), row=1, col=2)
    else:
        fig = go.Figure(data=go.Heatmap(
            z=matrix, x=x, y=list(range(n_var)),
            colorscale=colorscale_name, zmin=zmin, zmax=zmax,
            hovertemplate=('pos: %{x:,}<br>%{customdata}<br>'
                           f'{metric}: %{{z:.3f}}<extra></extra>'),
            customdata=np.array([[y_labels[r]] * data.analysis_length
                                 for r in range(n_var)]),
            colorbar=dict(title=metric),
        ))
        fig.update_xaxes(title_text=data.get_xlabel())

    # Highlight a specific variant row
    if highlight_variant_idx is not None:
        # Find row in sorted order
        for ri, (_, row) in enumerate(sorted_df.iterrows()):
            if row.name == highlight_variant_idx:
                fig.add_hline(y=ri, line=dict(color='red', width=2),
                              opacity=0.7)
                # Also add vline at variant position
                vpos = row['display_pos']
                if vpos and vpos != 0:
                    fig.add_vline(x=vpos,
                                  line=dict(color='red', width=1.5, dash='dot'),
                                  opacity=0.5)
                break

    if n_var <= 80:
        fig.update_yaxes(
            tickvals=list(range(n_var)), ticktext=y_labels,
            tickfont_size=max(5, 10 - n_var // 15),
            autorange='reversed')
    else:
        fig.update_yaxes(
            title_text=f"Variants (n={n_var}, 5'→3')",
            autorange='reversed')

    dir_label = {'all': '', 'decrease': ' (decreasing only)',
                 'increase': ' (increasing only)'}
    fig.update_layout(
        height=max(400, n_var * 16),
        title_text=f"{short}: {metric}{dir_label.get(direction_filter, '')}",
        margin=dict(l=140 if n_var <= 80 else 80, r=20, t=60, b=40),
    )
    return fig


# ============================================================================
# Tab 4: MPRA Lollipop
# ============================================================================

def build_mpra_lollipop(data):
    if data.mpra_df is None or 'value' not in data.mpra_df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No MPRA data loaded (use --mpra)",
                           xref="paper", yref="paper", x=0.5, y=0.5,
                           showarrow=False, font=dict(size=16))
        return fig

    df = data.mpra_df.copy()
    is_sig = df.get('significant', pd.Series([True] * len(df)))
    if 'ref' in df.columns and 'alt' in df.columns:
        df['mutation'] = df['ref'] + '>' + df['alt']
    else:
        df['mutation'] = 'unknown'

    mutation_colors = {
        'A>C': '#1f77b4', 'A>G': '#ff7f0e', 'A>T': '#d62728',
        'C>A': '#2ca02c', 'C>G': '#e377c2', 'C>T': '#9467bd',
        'G>A': '#8c564b', 'G>C': '#bcbd22', 'G>T': '#17becf',
        'T>A': '#7f7f7f', 'T>C': '#aec7e8', 'T>G': '#ffbb78',
    }
    fig = go.Figure()
    for mut_type in sorted(df['mutation'].unique()):
        mask = df['mutation'] == mut_type
        sub = df[mask]
        color = mutation_colors.get(mut_type, '#999999')
        for _, row in sub.iterrows():
            sc = '#2ca02c' if is_sig.loc[row.name] else '#cccccc'
            fig.add_trace(go.Scatter(
                x=[row['position'], row['position']], y=[0, row['value']],
                mode='lines', line=dict(color=sc, width=0.8),
                showlegend=False, hoverinfo='skip',
            ))
        opacities = np.where(is_sig[mask].values, 1.0, 0.3)
        fig.add_trace(go.Scatter(
            x=sub['position'], y=sub['value'], mode='markers',
            marker=dict(color=color, size=6,
                        line=dict(width=0.5, color='white'),
                        opacity=opacities),
            name=mut_type,
            hovertemplate=f'pos: %{{x:,}}<br>log₂: %{{y:.2f}}<br>{mut_type}<extra></extra>',
        ))

    fig.add_hline(y=0, line=dict(color='black', width=0.5))
    fig.update_layout(
        height=400, title_text='MPRA Variant Effect (Log₂)',
        xaxis_title=data.get_xlabel(), yaxis_title='Log₂ variant effect',
        legend=dict(orientation='h', yanchor='bottom', y=1.02),
        margin=dict(l=60, r=20, t=60, b=40),
    )
    return fig


# ============================================================================
# Tab 5: Comparison (with vlines and standardized axes)
# ============================================================================

def build_comparison(data, idx1, idx2, bin_label):
    x = data.get_display_positions()
    short = data.BIN_SHORT.get(bin_label, bin_label)
    wt_occ = data.reorient(data.wt_occ.get(bin_label,
                            np.zeros(data.analysis_length)))

    fig = make_subplots(
        rows=2, cols=3, shared_xaxes=True,
        horizontal_spacing=0.07, vertical_spacing=0.08,
        column_widths=[0.38, 0.34, 0.28])

    colors = ['#e85d50', '#2ca02c']
    all_delta_max = []

    for row_i, vidx in enumerate([idx1, idx2]):
        r = row_i + 1
        row_info = data.summary.iloc[vidx]
        hkey = row_info['hdf5_key']
        vdata = data.get_variant_data(hkey)
        color = colors[row_i]
        label_text = row_info['display_label']
        vpos = row_info['display_pos']

        if vdata is None or bin_label not in vdata:
            continue

        ld = vdata[bin_label]
        all_delta_max.append(np.max(np.abs(ld['delta_obs'])))

        fig.add_trace(go.Scatter(
            x=x, y=wt_occ, mode='lines',
            line=dict(color='#4a86c8', width=1.2),
            name='WT', showlegend=(row_i == 0), legendgroup='wt',
        ), row=r, col=1)
        fig.add_trace(go.Scatter(
            x=x, y=ld['variant_occ'], mode='lines',
            line=dict(color=color, width=1.2),
            name=label_text, showlegend=True,
        ), row=r, col=1)

        fig.add_trace(go.Scatter(
            x=x, y=ld['delta_obs'], mode='lines',
            line=dict(color=color, width=1),
            fill='tozeroy', fillcolor=_hex_to_rgba(color, 0.15),
            showlegend=False,
        ), row=r, col=2)
        fig.add_hline(y=0, line=dict(color='gray', width=0.5), row=r, col=2)

        nlp = -np.log10(np.maximum(ld['empirical_p'], 1e-10))
        fig.add_trace(go.Scatter(
            x=x, y=nlp, mode='lines',
            line=dict(color=color, width=0.8), showlegend=False,
        ), row=r, col=3)
        fig.add_hline(y=-np.log10(0.05),
                      line=dict(color='orange', width=0.5, dash='dash'),
                      row=r, col=3)

        # Variant position dotted line
        if vpos and vpos != 0:
            for c in [1, 2, 3]:
                fig.add_vline(x=vpos,
                              line=dict(color=color, width=1.5, dash='dot'),
                              opacity=0.5, row=r, col=c)

    if all_delta_max:
        dmax = max(all_delta_max) * 1.1
        for r in [1, 2]:
            fig.update_yaxes(range=[-dmax, dmax], row=r, col=2)

    for c in [1, 2, 3]:
        fig.update_xaxes(title_text=data.get_xlabel(), row=2, col=c)
    fig.update_layout(
        height=500, title_text=f'{short}: Variant Comparison',
        margin=dict(l=50, r=20, t=60, b=40),
    )
    return fig


# ============================================================================
# Layout
# ============================================================================

def build_layout(data):
    n_var = len(data.summary)

    sort_options = [
        {'label': "Position (5' → 3')", 'value': 'pos_asc'},
        {'label': "Position (3' → 5')", 'value': 'pos_desc'},
        {'label': 'Significance (best p)', 'value': 'sig'},
        {'label': 'Effect size (max |Δ|)', 'value': 'effect'},
        {'label': 'Coverage (high → low)', 'value': 'coverage'},
        {'label': 'FDR q-value', 'value': 'fdr'},
    ]

    sorted_summary = data.summary.sort_values('display_pos')
    variant_options = [
        {'label': f"{row['display_label']} (n={row['n_reads']}, q={row['fdr_q']:.3f})",
         'value': int(orig_idx)}
        for orig_idx, row in sorted_summary.iterrows()
    ]

    bin_options = [{'label': data.BIN_SHORT.get(l, l), 'value': l}
                   for l in data.bin_labels]
    heatmap_metrics = [
        {'label': '-log₁₀(p)', 'value': '-log10(p)'},
        {'label': 'Δ (signed)', 'value': 'Δ (signed)'},
        {'label': '|Δ|', 'value': '|Δ|'},
        {'label': 'Z-score', 'value': 'Z-score'},
    ]
    delta_modes = [
        {'label': 'Gray non-significant', 'value': 'gray'},
        {'label': 'Blank non-significant', 'value': 'blank'},
        {'label': 'Show all', 'value': 'all'},
    ]
    direction_options = [
        {'label': 'All changes', 'value': 'all'},
        {'label': 'Decreasing only', 'value': 'decrease'},
        {'label': 'Increasing only', 'value': 'increase'},
    ]

    # ── Tab 1: WT Landscape ───────────────────────────────
    tab1 = html.Div([
        dcc.Checklist(
            id='wt-show-mpra',
            options=[{'label': ' Show MPRA overlay', 'value': 'mpra'}],
            value=['mpra'] if data.mpra_df is not None else [],
            inline=True, style={'marginBottom': '5px'}),
        dcc.Graph(id='wt-landscape'),
    ], style={'padding': '10px'})

    # ── Tab 2: Variant Browser ────────────────────────────
    tab2 = html.Div([
        html.Div([
            html.Div([
                html.Label('Sort by:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='sort-dropdown', options=sort_options,
                             value='pos_asc', clearable=False, style={'width': '220px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center', 'marginRight': '20px'}),
            html.Div([
                html.Label('Δ display:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='delta-mode', options=delta_modes,
                             value='gray', clearable=False, style={'width': '200px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center'}),
        ], style={'marginBottom': '10px'}),
        html.Div([
            html.Button('◀ Prev', id='prev-btn', n_clicks=0,
                        style={'marginRight': '10px', 'padding': '5px 12px'}),
            dcc.Slider(id='variant-slider', min=0, max=n_var - 1,
                       step=1, value=0, marks=None,
                       tooltip={'placement': 'bottom', 'always_visible': False},
                       updatemode='mouseup'),
            html.Button('Next ▶', id='next-btn', n_clicks=0,
                        style={'marginLeft': '10px', 'padding': '5px 12px'}),
        ], style={'display': 'flex', 'alignItems': 'center', 'marginBottom': '5px'}),
        html.Div(id='variant-info',
                 style={'fontSize': '14px', 'fontWeight': 'bold',
                        'marginBottom': '5px', 'padding': '8px',
                        'backgroundColor': '#f8f9fa', 'borderRadius': '4px'}),
        dcc.Graph(id='variant-plot'),
        dcc.Markdown(id='cluster-details',
                     style={'padding': '10px', 'backgroundColor': '#f8f9fa',
                            'borderRadius': '4px', 'fontSize': '13px'}),
        dcc.Store(id='sorted-indices', data=list(range(n_var))),
    ], style={'padding': '10px'})

    # ── Tab 3: Library Heatmap ────────────────────────────
    tab3 = html.Div([
        html.Div([
            html.Div([
                html.Label('Bin:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='heatmap-bin', options=bin_options,
                             value=data.bin_labels[0], clearable=False, style={'width': '150px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center', 'marginRight': '15px'}),
            html.Div([
                html.Label('Metric:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='heatmap-metric', options=heatmap_metrics,
                             value='-log10(p)', clearable=False, style={'width': '150px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center', 'marginRight': '15px'}),
            html.Div([
                html.Label('Direction:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='heatmap-direction', options=direction_options,
                             value='all', clearable=False, style={'width': '170px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center', 'marginRight': '15px'}),
            html.Div([
                html.Label('Colors:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='heatmap-colorscale', options=COLORSCALE_OPTIONS,
                             value=None, clearable=True, placeholder='Auto',
                             style={'width': '170px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center'}),
        ], style={'marginBottom': '8px', 'display': 'flex', 'flexWrap': 'wrap'}),

        html.Div([
            html.Label('Scale range:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
            dcc.RangeSlider(id='heatmap-scale', min=-5, max=5, step=0.1,
                            value=[-5, 5],
                            marks={-5: '-5', -2: '-2', 0: '0', 2: '2', 5: '5'},
                            tooltip={'placement': 'bottom', 'always_visible': False}),
        ], style={'marginBottom': '8px', 'maxWidth': '600px'}),

        dcc.Checklist(
            id='heatmap-use-custom-scale',
            options=[{'label': ' Use custom scale range', 'value': 'custom'}],
            value=[], inline=True, style={'marginBottom': '8px'}),

        dcc.Graph(id='heatmap-plot'),
    ], style={'padding': '10px'})

    # ── Tab 4: MPRA ───────────────────────────────────────
    tab4 = html.Div([
        dcc.Graph(id='mpra-lollipop', figure=build_mpra_lollipop(data)),
    ], style={'padding': '10px'})

    # ── Tab 5: Comparison ─────────────────────────────────
    tab5 = html.Div([
        html.Div([
            html.Div([
                html.Label('Variant 1:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='cmp-var1', options=variant_options,
                             value=variant_options[0]['value'], clearable=False,
                             style={'width': '350px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center', 'marginRight': '15px'}),
            html.Div([
                html.Label('Variant 2:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='cmp-var2', options=variant_options,
                             value=variant_options[min(1, len(variant_options)-1)]['value'],
                             clearable=False, style={'width': '350px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center', 'marginRight': '15px'}),
            html.Div([
                html.Label('Bin:', style={'fontWeight': 'bold', 'marginRight': '8px'}),
                dcc.Dropdown(id='cmp-bin', options=bin_options,
                             value=data.bin_labels[0], clearable=False,
                             style={'width': '150px'}),
            ], style={'display': 'inline-flex', 'alignItems': 'center'}),
        ], style={'marginBottom': '10px', 'display': 'flex', 'flexWrap': 'wrap'}),
        dcc.Graph(id='comparison-plot'),
    ], style={'padding': '10px'})

    return html.Div([
        html.H2('MPRA Fiber-seq Library Browser',
                style={'textAlign': 'center', 'padding': '10px', 'marginBottom': '0'}),
        html.P(f'{n_var} variants | {data.analysis_length} bp | {data.get_xlabel()}',
               style={'textAlign': 'center', 'color': 'gray', 'marginTop': '0'}),
        dcc.Tabs([
            dcc.Tab(label='WT Landscape', children=tab1),
            dcc.Tab(label='Variant Browser', children=tab2),
            dcc.Tab(label='Library Heatmap', children=tab3),
            dcc.Tab(label='MPRA Data', children=tab4),
            dcc.Tab(label='Comparison', children=tab5),
        ]),
    ])


# ============================================================================
# Callbacks
# ============================================================================

def register_callbacks(app, data):

    # ── WT Landscape ──────────────────────────────────────
    @app.callback(
        Output('wt-landscape', 'figure'),
        Input('wt-show-mpra', 'value'),
        Input('wt-landscape', 'hoverData'),
    )
    def update_wt(show_mpra_val, hover_data):
        hover_pos = None
        if hover_data and 'points' in hover_data:
            pts = hover_data['points']
            if pts:
                hover_pos = pts[0].get('x')
        return build_wt_landscape(data,
                                   show_mpra='mpra' in (show_mpra_val or []),
                                   hover_pos=hover_pos)

    # ── Sort ──────────────────────────────────────────────
    @app.callback(
        Output('sorted-indices', 'data'),
        Output('variant-slider', 'value'),
        Input('sort-dropdown', 'value'),
    )
    def update_sort(sort_key):
        df = data.summary
        if sort_key == 'pos_asc':
            order = df['display_pos'].argsort().tolist()
        elif sort_key == 'pos_desc':
            order = df['display_pos'].argsort()[::-1].tolist()
        elif sort_key == 'sig':
            order = df['best_cluster_p'].argsort().tolist()
        elif sort_key == 'effect':
            max_d = np.zeros(len(df))
            for label in data.bin_labels:
                short = data.BIN_SHORT.get(label, label)
                col = f'{short}_max_delta'
                if col in df.columns:
                    max_d = np.maximum(max_d, df[col].values)
            order = (-max_d).argsort().tolist()
        elif sort_key == 'coverage':
            order = (-df['n_reads'].values).argsort().tolist()
        elif sort_key == 'fdr':
            order = df['fdr_q'].argsort().tolist()
        else:
            order = list(range(len(df)))
        return order, 0

    # ── Prev/Next ─────────────────────────────────────────
    @app.callback(
        Output('variant-slider', 'value', allow_duplicate=True),
        Input('prev-btn', 'n_clicks'),
        State('variant-slider', 'value'),
        prevent_initial_call=True,
    )
    def prev_variant(n, current):
        return max(0, current - 1)

    @app.callback(
        Output('variant-slider', 'value', allow_duplicate=True),
        Input('next-btn', 'n_clicks'),
        State('variant-slider', 'value'),
        prevent_initial_call=True,
    )
    def next_variant(n, current):
        return min(len(data.summary) - 1, current + 1)

    # ── Variant plot ──────────────────────────────────────
    @app.callback(
        Output('variant-plot', 'figure'),
        Output('cluster-details', 'children'),
        Output('variant-info', 'children'),
        Input('variant-slider', 'value'),
        Input('delta-mode', 'value'),
        State('sorted-indices', 'data'),
    )
    def update_variant_plot(slider_val, delta_mode, sorted_indices):
        if sorted_indices is None or slider_val >= len(sorted_indices):
            return go.Figure(), "", ""
        actual_idx = sorted_indices[slider_val]
        row = data.summary.iloc[actual_idx]
        fig, cluster_md = build_variant_plot(data, actual_idx, delta_mode)
        info = (f"#{slider_val + 1}/{len(sorted_indices)}  |  "
                f"{row['display_label']}  |  n={row['n_reads']}  |  "
                f"p={row['best_cluster_p']:.4f}  |  FDR q={row['fdr_q']:.4f}")
        return fig, cluster_md, info

    # ── Heatmap ───────────────────────────────────────────
    @app.callback(
        Output('heatmap-plot', 'figure'),
        Input('heatmap-bin', 'value'),
        Input('heatmap-metric', 'value'),
        Input('heatmap-direction', 'value'),
        Input('heatmap-colorscale', 'value'),
        Input('heatmap-scale', 'value'),
        Input('heatmap-use-custom-scale', 'value'),
        Input('mpra-lollipop', 'clickData'),
    )
    def update_heatmap(bin_label, metric, direction, colorscale,
                       scale_range, use_custom, mpra_click):
        zmin = scale_range[0] if 'custom' in (use_custom or []) else None
        zmax = scale_range[1] if 'custom' in (use_custom or []) else None

        # Find highlighted variant from MPRA click
        highlight_idx = None
        if mpra_click and 'points' in mpra_click:
            pts = mpra_click['points']
            if pts:
                click_pos = pts[0].get('x')
                if click_pos is not None:
                    # Find variant closest to this position
                    diffs = np.abs(data.summary['display_pos'].values - click_pos)
                    closest = np.argmin(diffs)
                    if diffs[closest] < 3:
                        highlight_idx = closest

        return build_heatmap(data, bin_label, metric,
                             direction_filter=direction,
                             zmin_override=zmin, zmax_override=zmax,
                             colorscale_name=colorscale,
                             highlight_variant_idx=highlight_idx)

    # ── Comparison ────────────────────────────────────────
    @app.callback(
        Output('comparison-plot', 'figure'),
        Input('cmp-var1', 'value'),
        Input('cmp-var2', 'value'),
        Input('cmp-bin', 'value'),
    )
    def update_comparison(idx1, idx2, bin_label):
        return build_comparison(data, idx1, idx2, bin_label)


# ============================================================================
# Main
# ============================================================================

def main():
    args = parse_args()
    print("=" * 60)
    print("MPRA Fiber-seq Library Browser")
    print("=" * 60)
    print(f"Loading: {args.h5}")

    data = LibraryData(
        args.h5, mpra_path=args.mpra,
        chrom=args.chrom, genomic_start=args.genomic_start,
        genomic_end=args.genomic_end, strand=args.strand,
    )

    app = Dash(__name__)
    app.layout = build_layout(data)
    register_callbacks(app, data)

    print(f"\nStarting server at http://{args.host}:{args.port}")
    print("Press Ctrl+C to stop.\n")
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == '__main__':
    main()
