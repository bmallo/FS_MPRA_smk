# Truth-set figure — analysis log

## 2026-11-05 — definitions, code, test on census elements

- User chose: confident ClinVar = >= 2 stars; arm B universe = all possible SNVs; likely pathogenic = rare (absent or joint AF < 0.01%) + PromoterAI |score| >= 0.2 + AlphaGenome large change; elements = all except GATA2 (replaced by TERT). Written to `rules.json` (with the immune and TERT AlphaGenome track sets) before any immune/TERT AlphaGenome score was seen.
- `evidence.py` (per-SNV evidence, counts, funnel) and `figure.py` (panels a–c). Tested on LDLR, HBB, SLC39A4 only (scratch output, not results/): likely benign candidates 0 / 0 / 3, matching the census truth sets after removing the ClinVar >= 2-star B/LB SNVs (LDLR rs17249134 and HBB rs139703273 move to arm A).
- Waiting for AlphaGenome (TERT, then the five immune promoters; started 2026-11-05 ~12:00, ~1.1 s per variant).

## 2026-11-05 — draft 1 without AlphaGenome

- User asked for a first draft without AlphaGenome; rule `draft1_no_alphagenome` written first. `evidence.py --no-alphagenome` and `figure.py --no-alphagenome` -> `results/truth_set_figure_20261105/draft1_no_alphagenome/`.
- Bug found and fixed before the draft was used: ClinVar keys for the immune windows took the chromosome of the first window (chr7) for all five, so JAK3 and STAT1 confident B/LB were missed. Now per window; all 48 confident ClinVar SNV records match an SNV in the table.
- Arm A (>= 2 stars) P/LP: LDLR 4, HBB 21, others 0. B/LB: LDLR 5, HBB 1, SLC39A4 1, JAK3 2, STAT1 1, TERT 13.
- Arm B without AlphaGenome: likely pathogenic 22–383 per element (1,771 pooled); likely benign 0–9 (24 pooled; 22 never tested in GTEx).

## 2026-11-05 — draft 2 with AlphaGenome

- AlphaGenome finished: TERT 913 and immune 15,145 variant-contexts, 0 failed requests, 0 duplicate rows; 0 SNVs missing a score in the evidence table. Track sets as rules.json (immune 7 biosamples, TERT 5, census as before).
- `evidence.py` + `figure.py` (full rules) -> `results/truth_set_figure_20261105/draft2_alphagenome/`. Arm A unchanged. Likely pathogenic 22–293 per element (1,368 pooled; 1,771 without AlphaGenome). Likely benign 0–8 (22 pooled; 24 without AlphaGenome; HBB 1 -> 0, CARD11 9 -> 8).

## 2026-11-05 — draft 3: >= 1-star ClinVar, benign = in gnomAD at any frequency, no GTEx

- User asked: ClinVar down to 1 star; benign needs only presence in gnomAD; drop GTEx. Rule `draft3` written first. `evidence.py --min-stars 1 --benign in_gnomad` (new options; defaults reproduce draft 2's counts exactly) and `figure.py` (reads the settings from summary.json; GTEx row removed from panel b) -> `results/truth_set_figure_20261105/draft3/`.
- Arm A P/LP SNVs: LDLR 8, HBB 24 (others 0); B/LB: LDLR 8, HBB 3, SLC39A4 2, JAK3 2, STAT1 4, TERT 22. Confident indels 7.
- Likely benign 72–389 per element (1,654 pooled, was 22); likely pathogenic 1,365 pooled (was 1,368; 3 SNVs moved into arm A).
