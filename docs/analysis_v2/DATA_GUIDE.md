# Working with the analysis v2 data

How to load the data, rerun the tests, and rebuild the explorer and the deck. Every recipe here was run on
2026-11-02 and reproduced the stored results exactly.

## Setup

- Compute on the standing node (`klone-node`; n3446 at the time of writing, 4 CPUs, 200 GB) or in Slurm,
  never on a login node.
- Run from the repo root:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=.
PY=/mmfs1/home/bmallo/miniconda3/envs/general/bin/python
```

- Tests: `$PY -m pytest -q tests/test_analysis_v2_*.py tests/test_explorer_data.py` (68 pass, ~15 s).
- Scripts in `tools/analysis_v2/` load each other with `importlib` (no package `__init__`); `from
  tools.analysis_v2 import access` works with `PYTHONPATH=.`.

## Coordinates

| System | HepG2 (contig `LDLR_MPRA_NotI`) | HEK293T (contig `LDLR`) |
|---|---|---|
| Element on the plasmid (0-based, half-open) | 1234–1445, forward | 3236–3447, **reversed** |
| Oligo (Kircher 2019; SNV names) | `pos - 1234 + 21` | `3446 - pos + 21` |

The element spans oligo 21–231 and has 101 A/Ts (one methylatable adenine each). `access.oligo_of()` does
the conversion. `footprints.tsv` columns `mid_*` are element indices (0-based); `oligo_*` and `span_*` are oligo.

## Inputs outside `results/analysis_v2_20261028`

| What | Path |
|---|---|
| Run registry (BAMs, contigs, reaction) | `config/analysis_v2/runs.json` (paths relative to `/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/fiberseq`) |
| Every setting, result and decision, by step | `config/analysis_v2/rules.json` |
| HepG2 run-3 ledger and design catalog | `results/rep3_ledger_20261016/`, `.../catalog/catalog.json` |
| LLR emissions (fold A) | HepG2 `results/llr_vs_fiberhmm_20261026/results/emissions.tsv.gz`; HEK `results/hek293t_llr_20261028/results/emissions.tsv.gz` |
| Published MPRA activity (Kircher 2019) | `/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/mpra_published/GRCh38_LDLR.tsv`, `GRCh38_LDLR.2.csv` |
| JASPAR PWMs | `extras/data/JASPAR2024_CORE_vertebrates_nr.jaspar` |

## Results layout (`results/analysis_v2_20261028/`)

Never overwrite or delete these; write new outputs to a new directory.

| Step | Directory | Key files | Tool |
|---|---|---|---|
| A1 | `a1_inventory/` | `inventory.tsv`, `differences.tsv` | `a1_inventory.py` |
| A2 | `a2_accounting/`, `a2_ledgers/` | `*/funnel.tsv`, `clone_depth_by_run.tsv` | `a2_*.py` |
| A3 | `a3_genotypes/` | discordance, collisions, second genotyper | `a3_genotype_audit.py` |
| A4 | `a4_informative/` | `visibility/variants_visibility.tsv` (per-variant annotations), `HepG2/positions.tsv` | `a4_*.py` |
| A5 | `a5_m6a_qc/` | `summary/` | `a5_*.py` |
| A6 | `a6_proxies/` | `extract/<run>/molecules.npz` (WT sample), `summary/` | `a6_*.py` |
| A7 | `a7_snv_edits/` | per-run edit effects | `a7_snv_edits.py` |
| A8 | `a8_confidence/` | `spikes/`, `repro/` | `a8_*.py` |
| A9 | `a9_statistics/` | **`extract/<run>/`** (per-molecule table for all later steps), `<run>/mde.parquet`, `calibration.tsv`, `balance.tsv` | `a9_extract.py`, `a9_framework.py` |
| A10 | `a10_effects/` | `HepG2_run3{,_pm5}{,_empirical_null}/` `effects.parquet`, `hits.tsv`, `snvs.tsv`; `HEK293T_run1/motifs.tsv` | `a10_effects.py`, `a10_sensitivity.py`, `a10_followup.py` |
| A11 | `a11_motifs_function/` | `function_snvs.tsv`, `function.tsv`, `dose_snvs.tsv`, `dose.tsv` | `a11_motifs_function.py` |
| 60–86 | `fp60_86/` | `adenines.tsv`, `labeling.tsv`, `variants.tsv` | `fp60_86.py` |
| 114–138 | `fp114_138/` | same | `fp_region.py` |
| Footprints | `fp_features/` | `footprints.tsv`, `calls.npz`, `effects.parquet`, `wt_map.npz` | `fp_features.py` |
| Explorer data | `explorer_v6/` (current; v2–v5 older) | `base.json`, `chunks/hepg2_NN.json` | `explorer_data.py` |

Each step's figures are in its `figures/` subdirectory.

## The main tables

**Per-molecule extract** (`a9_statistics/extract/HepG2_run3`, `HEK293T_run1`): analysis-ready WT + SNV
molecules (A3/A6 read filters; HepG2 WT = all molecules of a 25% sample of WT clones).
- `molecules.npz`: `pos` (2,169 plasmid A/T positions), `observed` and `methylated` (bit-packed, molecules ×
  positions), `distal_frac` (labeling), `passes`, `strand`.
- `meta.npz`: `label` (`WT` or SNV name), `cls` (`WT`/`SNV`), `clone`, `name_crc`, `read_quality`.
- `clones.tsv`: clone, label, readable barcode, molecules. No read names are stored.

**A10 per-feature effects** (`a10_effects/HepG2_run3/effects.parquet`, 598 k rows): one row per SNV ×
endpoint × feature. `endpoint` E1 (feature = element adenine, `bin` = adenine) or E2 (feature = position ×
size bin). `D` = variant − WT (fraction, labeling-standardised), `se` (random effects), `p` (NaN = not
testable or masked), `z`, `n` molecules, `C` clones, `expected` minority events, `masked`, `mde`, `hit`
(hierarchical FDR), `oligo`, `snv_oligo`, `caller_class`, `group`. `hits.tsv` = hit rows plus `distance`.
`snvs.tsv` = per SNV and endpoint: Simes p, selected, hits, max |z|. The `_pm5` directory is the ±5 bp mask
run; `_empirical_null` rescales z by the robust SD (post hoc).

**Function table** (`a11_motifs_function/function_snvs.tsv`): per SNV `act`, `act2` (Kircher log2 effect,
LDLR and LDLR.2), `C1` element-wide shift and `C2` local 11–50 bp change (points of inaccessible fraction,
±5 run), `rms`, `far_hits`, `region` (the 17 Repeat 3 SNVs).

**Footprints** (`fp_features/`): `footprints.tsv` (F01–F24: box `oligo_lo..oligo_hi` × `len_lo..len_hi`,
`span_lo..span_hi`, `wt_occupancy`, `wt_occupancy_confident`, `confident_share`, occupancy by labeling,
`flag`); `calls.npz` (`rows` = indices into the extract for molecules with labeling ≥ 0.2; per molecule up to
32 calls: `mid`, `length`, `start`, `end` in oligo, −1 padding in `mid`; `S` = Σ log(1 − p_a) over the call's
unmarked adenines); `effects.parquet` (SNV × mask none/A7/pm5 × calls all/confident × footprint: `occ`, `D`,
`se`, `p`, `sig`).

## Recipes

### Load molecules and rerun the variant-vs-WT test (`tools/analysis_v2/access.py`)

```python
from tools.analysis_v2 import access
d = access.load_element('HepG2')        # ~7 s; d.X molecules x 101 adenines: -1 not observed, 0 inaccessible, 1 accessible
ref = access.E1Reference(d)             # WT per labeling stratum + WT clone variance
masks = access.a7_masks('HepG2')        # {snv: masked oligo positions}; pm5=True for the +/-5 bp mask
t = ref.test(d, '79C>A', mask=masks['79C>A'])   # oligo, n, C, D, se, p, z, testable, masked
wt = ref.wt_fraction()                  # WT inaccessible fraction per adenine
```

`ref.test` gives the same D, se and testability as `effects.parquet` (`tests/test_analysis_v2_access.py`).
For a new grouping (pool of SNVs, subset of molecules), build `Y = X == 0` and `O = (X >= 0) & (labeling > 0.01)`
for those molecules and call `a9_framework.contrast` as `E1Reference.test` does. HEK293T loads the same way,
but per-SNV tests are not valid there (one clone per SNV); use clone means (`a10_effects.clone_mean_test`).

### Footprint occupancy per molecule

```python
import numpy as np, pandas as pd
R = 'results/analysis_v2_20261028/'
z = np.load(R + 'fp_features/calls.npz')
meta = np.load(R + 'a9_statistics/extract/HepG2_run3/meta.npz')
L = np.load(R + 'a9_statistics/extract/HepG2_run3/molecules.npz')['distal_frac'][z['rows']]
fp = pd.read_csv(R + 'fp_features/footprints.tsv', sep='\t').set_index('feature')
conf = (np.log(1 - np.clip(L, .01, .99)) / np.log(1 - 0.70))[:, None] * z['S'] / np.log(10) < -3
b = fp.loc['F04']
inbox = (z['mid'] >= b.oligo_lo) & (z['mid'] <= b.oligo_hi) & (z['length'] >= b.len_lo) & (z['length'] <= b.len_hi)
wt = meta['cls'][z['rows']] == 'WT'
inbox.any(1)[wt].mean(), (inbox & conf).any(1)[wt].mean()   # 0.0261, 0.0225 = footprints.tsv
```

These calls are unmasked. For masked calls on a variant, set the masked adenines to −1 in both the variant's
and WT molecules and call again with `fp_features.call_table` (see `fp_features._snv`).

### Common reads of the stored results

```python
eff = pd.read_parquet(R + 'a10_effects/HepG2_run3/effects.parquet',
                      filters=[('endpoint', '==', 'E1'), ('snv', '==', '146C>T')])
hits = pd.read_csv(R + 'a10_effects/HepG2_run3/hits.tsv', sep='\t')   # all hits with distance from the SNV
fe = pd.read_parquet(R + 'fp_features/effects.parquet').query('mask == "A7" and calls == "confident" and sig')
fn = pd.read_csv(R + 'a11_motifs_function/function_snvs.tsv', sep='\t')  # activity vs chromatin per SNV
```

## Conventions for new analyses

- Settings go into `config/analysis_v2/rules.json` (a new section) before results are looked at; label
  anything chosen afterwards as post hoc there and on slides.
- Compare variant with WT within a run, labeling-standardised, with the same mask on both sides.
- Keep HepG2 and HEK293T separate in every table and figure (HEK footprints are K01–K23, HepG2 F01–F24).
- Say accessible (marked) / inaccessible (unmarked); name metrics by what they measure.
- Report numbers, not conclusions; report a null with its MDE.
- Don't edit existing code, tests or results to make something pass; don't commit `results/` or read IDs.

## Promoter explorer

- Page: `tools/analysis_v2/explorer.html` (Plotly from cdnjs, loads `base.json`, then `chunks/hepg2_NN.json`
  on demand). Published at https://claude.ai/artifact/JTJkTUZ5NtDK6WSVt4BYsH (Version 7).
- Data (`explorer_data.py` docstring has the format; fractions are integers × 10⁴):

```bash
$PY tools/analysis_v2/explorer_data.py --root results/analysis_v2_20261028 \
    --catalog results/rep3_ledger_20261016/catalog/catalog.json \
    --features results/analysis_v2_20261028/fp_features --procs 4 --output results/analysis_v2_20261028/explorer_v7
```

- Look at it before publishing (headless Firefox on the node; static frame, no hover):

```bash
$PY tools/analysis_v2/render_page.py --page tools/analysis_v2/explorer.html \
    --data results/analysis_v2_20261028/explorer_v6 --work <scratch>/render --out <scratch>/shot.png \
    [--variant 146C>T] [--hash hek293t] [--dark] [--width 1400 --height 5200]
```

- Republish with the Artifact tool: `url` = the explorer link, `file_path` = `tools/analysis_v2/explorer.html`,
  `files` = `{"base.json": "<data dir>/base.json", "chunks/hepg2_00.json": "<data dir>/chunks/hepg2_00.json", ...}`
  (all 22 chunks). Read the artifact first if this session did not publish it.

## Deck

- Source: `docs/analysis_v2/deck/` (`deck.json` + `slides/*.html`, one section per step; copied from the
  session scratchpad on 2026-11-02). Published (Slides type) at https://claude.ai/artifact/NYSJ7cHat24e2SRj4SErwG.
- Slide images are `/_blob/<id>` assets stored in the artifact; the source PNGs are in each step's
  `results/.../figures/`. Upload new figures as assets of that artifact and reference the returned URL.
- Before publishing, read the remote deck (someone may have edited it), then publish `deck.json` with
  `root` = the deck folder and `files` mapping `slides/<id>.html`.
