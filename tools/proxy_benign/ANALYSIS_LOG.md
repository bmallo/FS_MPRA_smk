# Proxy benign column — analysis log

## 2026-11-06 — setup and runs

- Rules (`rules.json`) written from the user's choices before coordinates or data: in gnomAD (any AF); tested in GTEx V10 and no eQTL for any gene/tissue (untested does not count); no large AlphaGenome change (top-20%); ClinVar >= 1-star P/LP/B/LB excluded as benchmark; PromoterAI and primates not used.
- `elements.tsv`: 16 elements from the user (GATA2 = +9.5 enhancer, SHH = ZRS). AlphaGenome track sets from each gene's Fiber-seq lines present in AlphaGenome (WTC11, Jurkat E6-1, K562, HAP-1 DNase only; HBB uses K562). SHH interval re-centred (`ag_center`) so the 1 Mb window holds SHH and the ZRS; `ag_score.py` now reads `ag_center` and keeps WTC11 and HAP-1 tracks. Checked: SHH RNA_SEQ rows are present.
- Pulls into `results/proxy_benign_20261106/`: gnomAD (429 rate limit after 10 windows; the other 6 retried one at a time), GTEx batched (5,011 variants), ClinVar, reference.
- AlphaGenome: Slurm array 43007877 (`ag_array.sbatch`, cpu-g2, 16 tasks, at most 8 at once), ~1.35 s per variant per task.
- ClinVar counts from our windows differ from the user's table for several genes (e.g. APC P/LP 7 any / 6 at >= 1 star vs 3; MLH1 B/LB 26 / 21 vs 17; RBM8A P/LP 2 vs 0). The table's own P/LP and B/LB values are kept as given; the proxy column excludes our >= 1-star ClinVar SNVs.

## 2026-11-06 — scores and table

- Array 43007877: 16/16 COMPLETED (67–76 min each); 1–7 transient gRPC errors per task, all succeeded on retry.
- `proxy.py`: 0 SNVs missing AlphaGenome, 0 unmatched benchmark SNVs, 0 missing GTEx lookups. Per element (3,000 possible SNVs; SHH 3,429): in gnomAD 177–796; no large AlphaGenome change 1,935–2,238 (~66–75%, set by the top-20% rule on two readouts); GTEx class 0–3 (28 SNVs tested in GTEx across all elements, 22 of them eQTLs for some gene; 6 pass).
- Levels (% of possible SNVs): level 1 57–72%, level 2 3–18%, level 3 5 SNVs in total (TERT 1, PAH 1, HBG1 3). Bars total 66–79%.
- `table_figure.py` -> `table_proxy_benign.png/.pdf` (percent) and `_counts` version. Table values transcribed from the user's figure.

## 2026-11-06 — v2: AlphaGenome only, 350 bp window (no rerun)

- User: gnomAD and GTEx not helpful; AlphaGenome-only proxy benigns in a 350 bp window, only without an AlphaGenome rerun. Rule `v2_alphagenome_350bp` written first.
- MANE Select TSSs from the UCSC mane track (Ensembl REST too slow, stopped). A TSS-centred window ran 7-145 bp past 12 of the 14 promoter elements (no scores there), and APC and NIPBL have no MANE TSS inside. User chose (recorded in `window_revised` before counting): TSS-side 350 bp slid inward to stay inside the element; APC promoter 1B TSS (NM_001127511.3); NIPBL element 3' edge; GATA2/SHH centred.
- `proxy_350.py` -> `v2_350bp/windows_final.tsv`, `snv_350bp.tsv`, `counts.tsv`; uses the 1 kb top-20% cutoffs from `proxy.py`. Every window 1,050 SNVs; TSS inside every window except NIPBL. Proxy benign 394 (TERT) to 627 (NIPBL), 37.5-59.7%.
- `table_figure.py --v2` -> `v2_350bp/table_proxy_benign_350bp.png/.pdf`.

## 2026-11-06 — v2 cutoff comparison: top 20% of the 350 bp window vs of the 1 kb element

- Rule `cutoff_350_comparison` written first. `compare_cutoffs.py` -> `v2_350bp/cutoff_comparison.tsv`, `snv_350bp_cutoffs.tsv`, `counts_cut350.tsv`; `table_figure.py --v2 --cut350` -> `table_proxy_benign_350bp_cut350.png/.pdf`.
- Proxy benigns per gene: 1 kb cutoff 394–627 (7,744 pooled); 350 bp cutoff 670–791 (11,690 pooled). 7,690 called by both, 54 only with the 1 kb cutoff (MLH1 34, RBM8A 17, NIPBL 3), 4,000 only with the 350 bp cutoff. The 350 bp cutoffs are higher for RNA in every gene and for DNase in 13 of 16 (large effects concentrate near the TSS), so the 350 bp version is nearly a superset.
- Decision (user): keep the 1 kb cutoff. Final table: `v2_350bp/table_proxy_benign_350bp.png/.pdf`; the cut350 version is kept as a comparison only.
- User edits to the table: PAH row removed (15 genes) and HEL 92.1.7 removed from the Fiber-seq accessible column (HBB: HUDEP-2; ANKRD26: HEK293T, HAP1, Jurkat, WTC-11, HUDEP-2). Proxy counts unchanged: HBB still uses K562 (HUDEP-2 not in AlphaGenome), ANKRD26 its non-erythroid lines.
- User edits: RBM8A row replaced by PAH (15 genes); APC mechanism shown as LOE; footnote removed and grant-sized layout added (`table_figure.py --v2 --compact` -> `v2_350bp/table_proxy_benign_350bp_compact.png/.pdf`, 7.0 x 3.45 in, 300 dpi). The proxy rule now lives only in rules.json and this log, so the grant legend must state it.
