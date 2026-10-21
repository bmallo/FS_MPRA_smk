# Klone per-base profiles for the top-hit figures (run 2026-10-20)

Descriptive step from `docs/klone_agent_prompt_hit_figures.md`. It reads the existing site-methylation extraction; no BAM pass and no new statistics.

## Provenance

- Git commit: `5cbcebb` (branch `marginal-null-model`); `git diff fc4731a -- config/ldlr_site_methylation.json fsmpra/site_extract.py fsmpra/io.py fsmpra/processing.py` empty.
- Extraction: `results/rep3_site_methylation_20261021/extract` (job 42090780; audit `complete` True, 2,659,690 molecules).
- Slurm job 42094982 (`--account=stergachislab --partition=cpu-g2`, 1 CPU, 16G, threads = 1), node n3432, 2:50 elapsed, sacct MaxRSS 3,182,500 KB. Tests: `Ran 135 tests`, `OK (skipped=1)`. Log: `hit_profiles_42094982.log` (this directory).
- Command lines:
  - `python -m tools.variant_hit_figures profiles --extract results/rep3_site_methylation_20261021/extract --output results/rep3_site_methylation_20261021/hit_profiles` (0:21.96 wall, 3,155,796 KB max RSS)
  - `python -m tools.variant_hit_figures figures --profiles results/rep3_site_methylation_20261021/hit_profiles --output results/rep3_site_methylation_20261021/hit_figures` (0:43.53 wall, 347,292 KB max RSS)
- Files: `profiles.npz`, `rasters.npz`, `audit.json` from `hit_profiles/` (aggregate per-base fractions and 150 anonymous molecule rasters per allele; no read IDs). `hit_figures/` is not committed; it is on Klone at `results/rep3_site_methylation_20261021/hit_figures/`.

## Step-2 consistency readout (verbatim)

Profile columns average per-base fractions over the 211 bases; effects columns are from `docs/evidence/site_methylation_mask4_20261021/effects/effects.tsv.gz`. No allele differs by more than 1 pp (all listed pairs agree to 0.01 pp).

```
allele	molecules	U	M	C	profile_dMSP_pp	effects_dMSP_pp	profile_dNUC_pp	effects_dNUC_pp	finite_bases_msp_C_std
WT	1243891	93030	266758	884103	+0.00	+nan	+0.00	+nan	211
146C>G	2721	230	563	1928	-4.00	-4.00	+4.05	+4.05	211
146C>T	2831	212	595	2024	-4.28	-4.28	+4.15	+4.15	211
149C>G	3726	307	760	2659	-2.93	-2.93	+3.39	+3.39	211
98C>G	1885	142	437	1306	-3.45	-3.45	+2.93	+2.93	211
129C>G	2933	202	630	2101	-1.07	-1.07	+0.46	+0.46	211
133C>G	1960	137	439	1384	+1.63	+1.63	-1.43	-1.43	211
56G>C	1627	134	312	1181	+1.53	+1.53	-1.44	-1.44	211
56G>T	2240	153	477	1610	+0.10	+0.10	-0.12	-0.12	211
58G>T	2309	181	481	1647	-0.36	-0.36	+0.37	+0.37	211
117A>T	1859	167	405	1287	+4.45	+4.45	-3.74	-3.74	211
117A>C	1233	93	272	868	+4.65	+4.65	-1.50	-1.50	211
115A>C	1869	155	406	1308	+4.22	+4.22	-2.81	-2.81	211
191G>T	4935	345	1076	3514	+0.84	+0.84	+0.62	+0.62	211
73T>C	1735	116	396	1223	-0.83	-0.83	+0.14	+0.14	211
125A>T	1927	145	430	1352	+1.38	+1.38	+1.07	+1.07	211
alleles with finite msp_C_std: 631 of 634
raster alleles: 23 molecules per allele: [150]
{
 "schema_version": "fsmpra.hit_profiles.v1",
 "extract": "results/rep3_site_methylation_20261021/extract",
 "config": "config/ldlr_site_methylation_mask4.json",
 "strata_edges": [
  0.04260780289769173,
  0.07854209095239639,
  0.12012320011854172,
  0.18326488137245178,
  0.2643737196922302,
  0.34240245819091797,
  0.4111909568309784,
  0.4758726954460144,
  0.5385010242462158
 ],
 "strata_weights": [
  0.09951216091337774,
  0.09939905192042103,
  0.10021230557977973,
  0.10030053059428597,
  0.1001206872954848,
  0.0999317952772471,
  0.09992953309738797,
  0.10057425435724118,
  0.09937190576211143,
  0.10064777520266303
 ],
 "groups": "0 U, 1 M, 2.. C strata",
 "cohort": "eligible and footprint evaluated",
 "raster_per_allele": 150,
 "seed": 20261021,
 "states": [
  "unannotated",
  "accessible (MSP)",
  "nucleosome-sized 85-250 bp",
  "long protection >250 bp",
  "TF-sized footprint 10-60 bp"
 ],
 "interpretation": "Per-base fractions; C_std is standardized over WT distal-labeling strata. Rasters are random C-class molecules ordered by element MSP fraction. Descriptive; one preparation."
}
```
