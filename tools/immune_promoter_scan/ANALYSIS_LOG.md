# Immune promoter scan — analysis log

## 2026-11-05 — gnomAD, GTEx, PromoterAI for CARD11, JAK1, JAK3, STAT1, STAT3

- Windows, TSSs and transcripts from the user (`windows.tsv`); rules written before any retrieval (`rules.json`; census definitions). All five genes are on the minus strand; each user TSS = PromoterAI `tss_pos` + 1 for that transcript.
- Raw pulls: `fetch.py gnomad`, `promoterai`, `gtex --batched` with `--windows-tsv tools/immune_promoter_scan/windows.tsv` into `results/immune_promoter_scan_20261105/`; then `scan.py`. Checks: 0 GTEx lookups missing, 0 reference mismatches (gnomAD, PromoterAI).
- PromoterAI covers TSS +/- 500 bp only: not scored for 367 bp of CARD11, 124 bp of JAK1, 362 bp of STAT3 (reported, not counted as no effect).
- AlphaGenome and ClinVar not yet run.

## 2026-11-05 — ClinVar added (context only)

- Rule added to `rules.json` (`clinvar`) before the query; census classes and small-variant rule; parsed with `truth_sets.clinvar_table`. `scan.py` now writes `clinvar_variants.tsv` and ClinVar columns in `counts.tsv`, `variants.tsv`, `common.tsv`, `all_snvs.tsv`. 0 reference mismatches.
- 26 records, all SNVs in the windows: CARD11 2 (no classification: "evidence only" somatic splicing submissions from the Rogan lab on TCGA tumours), JAK1 0, JAK3 5, STAT1 10, STAT3 9. No P/LP anywhere. B/LB 6 SNVs (JAK3 2, STAT1 4), VUS 16, conflicting 2 (STAT3).
- JAK3 B/LB rs7254346 (+14) and rs199713123 (+24) are common but fail B2 (eQTL) and PromoterAI (-0.17, +0.15) respectively; the 4 STAT1 B/LB are < 1% in every group. So no ClinVar B/LB is a candidate benign SNV.

## 2026-11-05 — TERT and AlphaGenome

- TERT promoter chr5:1294953-1295254 (`windows_tert.tsv`, MANE ENST00000310581, TSS 1295068) through the same pulls and `scan.py --windows-tsv` into `results/tert_promoter_scan_20261105/`; 0 lookup/reference mismatches. ClinVar 82 records (P/LP SNV 1, indel 3; B/LB SNV 25).
- AlphaGenome (`ag_score.py --windows-tsv`, KEEP extended with immune and TERT biosamples; `census_lib.context()` maps new windows to their own interval) running for TERT then the immune windows.
