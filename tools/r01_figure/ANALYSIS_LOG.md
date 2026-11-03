# R01 preliminary-data figure (LDLR plasmid Fiber-seq MPRA): analysis log

Brief: `docs/subprojects/r01_ldlr_figure_prompt.md`. Builds on analysis v2 (`docs/analysis_v2/`); no re-derivation.
Outputs: `results/r01_figure_20261102/`. Code: `tools/r01_figure/`. Target: NIH R01 preliminary-data figure (not a
journal; general figure defaults apply).

## 2026-11-02 Step 1: re-orient (no new statistics)

- Code: `tools/r01_figure/reorient.py` (tests `tools/r01_figure/test_reorient.py`, 8 pass), run at git 78a4f18
  (uncommitted script) on n3446.
- Output: `results/r01_figure_20261102/step1_reorient/` (`regions.tsv`, `snv_footprint.tsv`, `leads.tsv`,
  `repressor_scan.tsv`, `motifs.tsv`, `summary.json` with md5 of every input).
- Inputs (md5 in summary.json): fp_features footprints/effects, A11 function_snvs/dose, A10 HEK motifs, HEK
  extract clones.tsv, A4 visibility, fp60_86 variants, Kircher 2019 LDLR and LDLR.2.
- Reproduced stored counts: 89 significant SNV x footprint pairs (A7 mask, confident calls), 79 at F04; HEK 197
  SNVs observed, 105 with >= 300 molecules, 119 with >= 3 clones. Activity parsed here equals A11 act/act2 (test).
- Definitions used only for description (no ranking yet): size class from the footprint's length box
  (TF = box overlapping 17-45 bp); "activator-like" = footprint D and published activity same sign,
  "repressor-like" = opposite; repressor scan = SNVs with activity > 0.5 log2 and P < 0.05, footprints within 20 bp.
- Descriptive check (exploratory, not prespecified): Spearman of F04 D (confident calls) vs activity across SNVs in
  66-82: A7 -0.10 (n 51, p 0.47); pm5 -0.19 (p 0.18); across 60-86 pm5 -0.31 (n 81, p 0.006).
- Noted: brief says 244 HEK SNVs observed; the analysis extract has 197 (post-filter). fp60_86 reports changes as
  WT - variant (positive = loss); fp_features D is variant - WT.
- Open questions for the user: HEK293T activity data? use the A11 element-wide shift as the repressor panel?
  clinical variant list from the paper figure; OK to download ClinVar/LOVD?

## 2026-11-02 User answers after step 1
- No HEK293T activity data: panel 3 compares chromatin only.
- Repressor panel: use a variant whose TF-sized footprint and activity move in opposite directions (68C>T is
  acceptable), showing the contradicting alleles.
- Clinical list: from the user's paper; ClinVar/LOVD download approved.

## 2026-11-02 Step 2: shortlist
- Rules written first: `tools/r01_figure/selection_rules.json` (before `rank.py` was run; step-1 tables had been
  viewed). Code `tools/r01_figure/rank.py`; output `results/r01_figure_20261102/step2_shortlist/ranked.tsv`.
- Attrition: 89 significant pairs (A7, confident) -> 87 TF-sized -> 66 pass the activity gate (P < 0.05 in LDLR,
  same sign in LDLR.2): 40 activator-like, 26 repressor-like.
- Explorer views (render_page.py, headless Firefox; absolute --out path needed, a relative one writes no file):
  79C>A, 71T>C, 71T>A, 68C>T, 149C>T, 98C>G in `step2_shortlist/explorer_views/`.

## 2026-11-02 Step 3 groundwork: clinical variants
- Paper: the Molecular Cell page returned 403; the variant list comes from the PMC/preprint version (PMC12154926),
  Fig. 6b + Ext. Data Fig. 12: c.-223A>G, c.-214T>C, c.-214T>A, c.-156C>T, c.-193C>T, c.-79C>G (the sixth Fig. 6b
  variant is unnamed there). Needs checking against the published version.
- Downloads (`results/r01_figure_20261102/clinical_raw/`, md5): ClinVar esearch e96422800e7360d952713b6eaa344886,
  esummary 1f910fb89c784f68f0b2f6a370de7492 (128 records); LOVD shared LDLR all aace3b3213e5a50f7e54c1a623ecb9ec
  (4,199 records; the window query returned [] so the full gene was taken and filtered).
- Mapping `tools/r01_figure/clinical.py` -> `step3_clinical/clinical_variants.tsv`: c. = oligo - 285, checked on
  40 ClinVar SNVs. 126 distinct variants in the element (111 SNV, 4 1-bp del, 2 1-bp ins, 9 larger); 117 designed.
  8 ClinVar records are outside or larger than the element (CNVs, c.-229_-90del, ...).

## 2026-11-02 Motif scan under F04 (user request; exploratory)
- Databases (latest releases, downloaded to `results/r01_figure_20261102/motif_db/`): JASPAR 2026 CORE vertebrates
  non-redundant (1,019 motifs, md5 c63218ac...), HOCOMOCO v14 H14CORE (1,595, md5 455a649b...; annotation jsonl
  183d473d...). Bundled JASPAR 2024 / HOCOMOCO v12 not used.
- Code `tools/r01_figure/motif_scan.py` (FIMO-style: 0.1 pseudocount, element background, exact P from the score
  distribution, checked against 400k simulated sequences: 1e-3 -> 9.7e-4 exact vs 9.3e-4 simulated).
- Output `results/r01_figure_20261102/motif_scan/`. WT element: 1,038,608 tests, 116 matches P < 1e-4, 0 with
  q < 0.05 (also true for the known SREBF1 and SP1 sites; the element is too short for element-wide q).
- F04 core 66-81 (GCCGATGTCACATCGG): 0 matches at P < 1e-4 (chance expectation 8.6); tiles 61-76 and 69-84 are the
  only zero tiles in the element (median 6). Positive controls found: SREBF1 127-136 P 2e-6; R3 Sp1 41 matches
  (SP4/SP1/KLF, P 2e-6); R1 8 matches (SP4, ZFP3).
- Alleles: 71T>C creates CGACGTCA (CRE-like): CREB1 P 7.3e-4 -> 3.1e-5, CREB3L4, ATF1, CREM. 71T>A strengthens SPIB/SPIC
  (1.3e-4 -> 6e-6). 68C>T: CREB1/ATF4 3.7e-4 -> 1.4e-5. Post hoc observations, not tests.

## 2026-11-02 Repeat 3 inventory (no new statistics)
- ClinVar P/LP SNVs in the element: 8. In R3 Sp1 (143-151): 149C>T, 150C>G (2); SRE-1/R2: 133C>T, 134C>G, 136C>A (3);
  other: 117A>G, 186A>G, 187C>T. Conflicting in R3: 145C>T, 146C>G. VUS at 138-157: 19.
- Per-SNV local Sp1-site footprint is mostly untestable: F15 (box 144-148, 17-23 bp, WT 0.25%) p = NaN for every
  SNV; A11 10-40 bp coverage testable only at 140-141, 146-147, 151, 153-155.

## 2026-11-02 Repeat 3 Sp1 site deep dive (approved plan; rules in `tools/r01_figure/repeat3_rules.json` first)
- Code `tools/r01_figure/repeat3.py` (tests `test_repeat3.py`; 14 pass in the folder), output
  `results/r01_figure_20261102/repeat3/` (105 s, 4 procs). Groups: G1 4, G2 21, G3 2, G4 33, G5 6, G6 6 SNVs.
- WT occupancy: site box (confident 10-40 bp calls, midpoint 143-151) 0.24% (all calls 0.57%); nucleosome-size box
  (confident >= 81 bp) 28.8%.
- Site box, pooled (pts, 95% CI): G1 -0.11 (-0.24, 0.01); G2 +0.05 (-0.10, 0.20); G3 +0.88 (0.49, 1.28); +/-5 mask
  ~0 for all. The local Sp1 footprint is at the detection floor: no measurable loss.
- Nucleosome-size, pooled: G1 +5.88 (3.32, 8.45); G2 +4.35 (3.60, 5.10); G3 +0.97 (-5.6, 7.5; k 2); G4 +0.51
  (-0.35, 1.37); G5 +0.14 (-0.82, 1.09); G6 -0.54 (-1.61, 0.53). Per SNV in the site: 26/27 > 0, 22 with z > 1.96.
  Within-site Spearman nucleosome change vs activity -0.69 (P 7e-5, n 27; exploratory). 147T>G (motif gain):
  -2.34 pts, activity +0.63.
- HEK293T check: access.load_element('HEK293T') WT inaccessible fraction matches explorer_v6 (max diff 5e-5).
  Clones >= 20 molecules: WT 147 clones (far residual 5-95%: -1.17 to +1.18 pts); site SNV clones 10, of which 7
  above the WT 95th percentile (+2.2 to +3.7); not above: 151C>T (+1.18), 147T>A (+0.58, activity-neutral in
  HepG2), 148C>A 21-molecule clone (-3.3).
- Draft figure `tools/r01_figure/fig_repeat3.py` -> `results/r01_figure_20261102/fig_repeat3_v3/` (pdf, png,
  panel_*.tsv). v1/v2 were layout drafts. Traced 149C>T nucleosome-size D independently from calls.npz:
  1,276 molecules, raw 31.19% vs WT 28.77%, standardised D 3.207 pts = stored 3.207.

## 2026-11-02 Review of the R01 Research Strategy draft (Specific_Aims_-_Plasmid_Fiber-seq_R01.docx, user upload)
- Read text and 25 comments (3 assigned to Ben: #20 analysis plan, #23 proxy-benign counts per element, #7 compact
  figure). Numbers checked against stored results: `results/r01_figure_20261102/grant_review/claim_numbers.txt`.
- Kircher LDLR file: 954 SNVs over oligo -33..284; HepG2 run 3 library covers oligo 21-231 = 633 SNVs; HEK extract
  197 SNVs. The draft's "950 of 954" does not match either v2 dataset.
- c.-156C>T (129C>T): A10 E2 gain in 20-40 bp coverage +1.4 to +2.0 pts over 120-142 and E1 +38 pts at 130, all
  within 13 bp, in the F11 region (3% confident calls); no fp_features footprint test available.
- ClinVar in element: 12 P/LP (8 SNV, 4 indel), 6 B/LB, 9 conflicting, 92 VUS. Sp1 site: 14 VUS among 25
  site-breaking SNVs; 13 of them nucleosome z > 1.96 and activity < -1 log2.

## 2026-11-02 gnomAD proxy-benign preview (grant Aim 1B; rules in `tools/r01_figure/gnomad_rules.json` first)
- gnomAD v4.1 GraphQL region download (`results/r01_figure_20261102/gnomad_raw/gnomad_r4_region.json`, md5
  0b2dc368977b91f7a0cc1ad8c11130ba): 227 SNV records, 120 with AC 0 dropped, 107 in the element (all designed),
  3 filtered. Only 68C>T has AF >= 1e-4 (AC 579, AF 1.1e-3, 3 homozygotes); none >= 1%.
- Code `tools/r01_figure/gnomad_check.py` (also computes the nucleosome-size endpoint for all 633 SNVs; 858 s).
  First run failed (worker pickling of a dynamically loaded module); fixed by registering it in sys.modules.
- Abnormal fractions (prespecified): activity P < 1e-5: P/LP 8/8, B/LB 3/6, gnomAD 50/107, not in gnomAD 235/526,
  all 285/633. Nucleosome |z| > 3.29: P/LP 1/8, B/LB 0/6, gnomAD 10/107, all 34/633. Footprint sig: P/LP 2/8,
  B/LB 2/6, gnomAD 14/107.
- POST HOC (not prespecified) effect-size view, `gnomad_check/posthoc_effect_size.txt`: P < 1e-5 and |log2| >= 1:
  P/LP 4/8, B/LB 1/6 (68C>T, a gain), gnomAD 24/107, not in gnomAD 101/526.
- Draft figure `tools/r01_figure/fig_gnomad.py` -> `results/r01_figure_20261102/fig_gnomad_v3/`.

## 2026-11-02 Predictor and database layer (LDLR only)
- GTEx (user wrote "GREx"; taken as GTEx): API v2. v8 ENSG00000130164.13: 74 significant single-tissue eQTL rows
  (7 tissues); v10 ENSG00000130164.14: 269 rows (13 tissues); no liver; none in the element (nearest ~1.5 kb).
  Element SNVs are presumably below GTEx's MAF >= 1% testing cutoff (max gnomAD AF 1.1e-3; inference, not checked
  against GTEx's variant list: the variantByLocation endpoint returned Not Found).
- PromoterAI: no public weights (repo has code only); precomputed scores need the academic license form (link by
  email). The repo's CAGI5_saturation.tsv has 777 LDLR rows with over/under/none labels but no scores.
- AlphaGenome: needs a user API key (deepmind.google.com/science/alphagenome); package not installed.
- Merged table `tools/r01_figure/merge_evidence.py` -> `results/r01_figure_20261102/evidence_v1/evidence.tsv`
  (633 SNVs x 30 columns; predictors pending). HEK293T: 188 SNVs with molecules, 114 with >= 3 clones.

## 2026-11-02 AlphaGenome scores and cross-evidence correlations (LDLR)
- AlphaGenome API (client 0.10.0, venv /mmfs1/gscratch/stergachislab/bmallo/large_home/envs/alphagenome; key in
  ~/.config/alphagenome/key, mode 600, never logged). Rules `alphagenome_rules.json` written first. Script
  `tools/r01_figure/alphagenome_score.py`: 633 SNVs, sequential, ~1 s each, 0 failures. Output
  `results/r01_figure_20261102/alphagenome/` (tracks.tsv md5 3bdbea46..., scores.tsv md5 23c7b371...).
- Merged `evidence_v2/evidence.tsv`. Correlation rules `correlation_rules.json` written before any AlphaGenome
  score was viewed; `tools/r01_figure/correlate.py` -> `correlations_v1/` (pairs.tsv, clinical.tsv).
- Prespecified (Spearman, all SNVs unless noted): MPRA vs element-wide shift -0.50 (n 631; Sp1 site -0.75);
  MPRA vs nucleosome-size -0.32 (site -0.69); MPRA vs local inaccessibility ~0; HepG2 vs HEK293T local 0.62 (n 113),
  element-wide/far 0.46; AlphaGenome RNA (HepG2) vs MPRA 0.12 (site 0.81, outside 0.05); AG CAGE vs MPRA 0.19;
  AG DNase vs nucleosome-size -0.45; AG SP ChIP vs nucleosome-size -0.42. AUROC P/LP (8) vs B/LB (6): MPRA 0.90,
  nucleosome 0.50, AG RNA 0.63, AG liver RNA 0.69, AG DNase 0.60.
- POST HOC (descriptive, `correlations_v1/spotlight.txt`): AG DNase vs MPRA 0.58. Region view: Sp1 R3 all
  readouts agree; SRE-1 MPRA only (median -0.95 log2; AG DNase ~0; Fiber-seq blind); F04 core: Fiber-seq footprint
  and MPRA alleles (71T>A/C), AG ~0; CRE 182-192: AG DNase -0.17, MPRA ~0.
- PromoterAI pending (user requested access).

## 2026-11-02 Figure option A (user choice)
- Display rules `tools/r01_figure/figA_rules.json` (panel b rule written before its fractions were computed; post hoc
  because region medians had been viewed). Script `tools/r01_figure/fig_A.py` -> `results/r01_figure_20261102/
  fig_A_v2/` (v1 = layout draft). Panel b: R3 Sp1 78/81/78% (MPRA/Fiber-seq/AlphaGenome); SRE-1 52/12/25%;
  F04 core 24/67/10%; CRE 0/6/64%; R1 56/23/48%; rest 10/12/11%.
- Traced (`fig_A_v2/trace.txt`): a 149C>T -2.62 = Kircher raw; b R3 MPRA 21/27 from raw file; c 149C>T 3.207 =
  calls.npz recompute; d rho 0.623, n 113 recomputed. Draft legend + caveats `fig_A_v2/legend_draft.md`.

## 2026-11-02 Panel a revision (user requests 1-4)
- TSS: UCSC ncbiRefSeq (`predictors_raw/refseq/ucsc_ncbiRefSeq_ldlr.json`, md5 deb4049c...): NM_000527.5 txStart
  chr19:11,089,463 = oligo 199 = c.-86; CDS start 11,089,549 = oligo 285 = c.1 (c.-1 = oligo 284).
- Conflicting variants: ClinVar VCV efetch (`clinical_raw/clinvar/conflicting_vcv.xml`, md5 f36ad371...). In-element
  SNVs: c.-140C>T 3 LP/1 VUS; c.-139C>G 5 P/LP/1 VUS; c.-121T>C 7 P/LP/5 VUS; c.-191C>A 2 LP/4 VUS; c.-97G>T
  1 VUS/1 LB; c.-84G>A 2 VUS/1 B.
- Display rules `tools/r01_figure/panel_a_rules.json` (post hoc display rule). Sites (`panel_a_sites/`): motif
  sites with P < 1e-5 and median MPRA <= -0.5: 89-100 (best ZFP3; SP4 only at 5e-5 -> "SP-like, weak"), 126-154
  (SREBF 127-136 and SP/KLF 143-151, drawn separately). Not drawn: ZNF77 23-46, TF7L1 167-178, CEBP/E2F 199-213
  (MPRA medians ~0). Motif-free confident TF-size footprints: F04 only (core 66-81 = c.-219..-204; 86% confident;
  0 motif matches at P < 1e-4; 13 ClinVar SNVs: 11 VUS, 1 B, 1 LB; 8 in gnomAD). F18 is confident (89%) but has
  16 motif matches.
- `fig_A.py` panel a rewritten (c. axis, TSS arrow, motif/footprint track, gnomAD rug, GTEx note) ->
  `results/r01_figure_20261102/fig_A_v4/` (panel_a.png/pdf standalone + full figure). Panel b row names not yet
  updated to the new labels.

## 2026-11-02 Panel a: gnomAD track dropped; confident-footprint version
- gnomAD rug made optional and turned off (not visually interpretable); gnomAD/GTEx go in the legend text.
- New `panel_a_footprints` (rule `panel_a_rules.json` footprint_version): confident (>= 0.8) TF-size footprints
  F04 (66-81 drawn, ~27 bp, 2.2% of molecules, no motif), F07 (84-108, ~23 bp, 0.25%, at SP-like site),
  F09 (84-130, ~45 bp, 0.081%), F18 (156-186, ~31 bp, 0.36%). Output `results/r01_figure_20261102/fig_A_v7/`.

## 2026-11-02 Panel a: motifs only (user decision)
- Motif-free footprint removed from panel a (showing one footprint next to motifs was confusing). Output
  `results/r01_figure_20261102/fig_A_v8/panel_a.png`. Panel b redesign (ClinVar variant mechanism categories) proposed,
  not yet run. Evidence available for the 104 ClinVar SNVs: nucleosome z 104, F04 test 102, any significant
  footprint 19, HEK >= 3 clones 25; local A/T (A4 window) 5-14.

## 2026-11-02 Panel b: mechanism categories for ClinVar SNVs (rules `classify_rules.json` written first)
- User addition: non-significant signed TF-size footprint differences count as supporting evidence (|mean z| >= 1).
- `tools/r01_figure/classify.py` -> `results/r01_figure_20261102/classify/` (variants.tsv for all 633, summary.tsv,
  notes.txt). First run failed (numpy array .fillna); fixed. 104 ClinVar SNVs: footprint evidence significant 37,
  supporting 20, none 47; breaks SP/KLF/SREBF motif 31; creates any motif 63 (permissive; decides 2 calls:
  166A>C, 210G>A); nucleosome q < 0.05 closes 21, opens 4; poorly visible 9.
- Primary (|log2| >= 0.5): P/LP 8 = 5 activator disrupted, 1 repressor lost (117A>G), 1 chromatin only (134C>G,
  -0.48 log2), 1 no effect (187C>T); conflicting 6 = 4 activator disrupted, 2 no effect; VUS 84 = 34 activator
  disrupted, 1 repressor gained (144T>C), 4 down unresolved, 4 activator gained, 16 chromatin only, 19 no effect,
  6 not assessable; B/LB 6 = 1 repressor lost (68C>T), 1 chromatin only (80G>T), 4 no effect.
- Activator disrupted (43) by region: Sp1 GC box 17, SRE 8, motif-free footprint 8, SP-like 6, other 4; by evidence:
  motif 27, significant footprint 11, supporting footprint 5. Mixed evidence (motif broken + significant local
  footprint gain): c.-156C>T, c.-153C>T, c.-152C>T (P/LP), c.-150A>G, c.-138T>C.
- Draft `tools/r01_figure/fig_panel_b.py` -> `results/r01_figure_20261102/panel_b_v2/`.

## 2026-11-02 c.-156C>T (129C>T) follow-up: created motifs and allele comparison (exploratory)
- Breakdown figure `tools/r01_figure/fig_variant.py` (v1, v2 in results/r01_figure_20261102/variant_c156CT_v*;
  v3 stopped by user request). Mask cache `results/r01_figure_20261102/cache/a7_masks.json` (complete, 633 x 2).
- `tools/r01_figure/motif_gain.py` -> `results/r01_figure_20261102/motif_gain_c156CT/`: 129C>T creates 6 motifs
  (variant P < 1e-4, >= 10-fold): ZN618 (P 1.2e-5), ZN584, ZN716, PO2F3, GLYR1, NFAT5 (+ PO2F2 7.5-fold); all overlap
  the protected stretch c.-165..-143. Base rate: 405/633 SNVs create >= 1 (median 1); 129C>T at the 91st percentile
  by count, 81st by best P.
- Allele comparison (`allele_comparison.txt`): at c.-156, C>T footprint gain z +4.7 (35 hits), C>A gain z +5.3
  (36 hits; creates 2 motifs, only GLYR1 shared), C>G loss z -2.3; expression -2.52 / -2.33 / -2.49. The gain tracks
  introducing an A/T at 129 (making a 4-bp A/T stretch c.-158..-155), not a specific created motif: consistent with a
  Hia5 sequence-context effect read as protection by WT-built emissions. Revises the earlier reading of the
  labeling plot (confident calls rising with labeling is also expected under that artifact).

## 2026-11-02 PromoterAI access
- User received the license; scores and models are in two Illumina OneDrive folders. Anonymous access from klone
  returns 403 and rclone has only S3 remotes, so files must be downloaded by the user and copied to
  /gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA/promoterai/ (created; group quota 93% used,
  ~15 TB free).
- User chose precomputed scores. `tools/r01_figure/promoterai_ingest.py` written (window filter, ref check, LDLR/+
  strand preference, tabix if indexed); passes a toy-file check. Waiting for the file on klone.

## 2026-11-02 PromoterAI precomputed scores
- Files (user copy from Illumina OneDrive): projects/LDLR_MPRA/promoterai/promoterAI_tss500.tsv.gz (md5
  f6ac61c2ea66313222c4ba2d1c94033e; plain gzip, columns chrom pos ref alt gene gene_id transcript_id strand tss_pos
  promoterAI) and tss.tsv (md5 f3fdbd4410b2573eb8b54dbf0622ffb6).
- `promoterai_ingest.py` (added awk prefilter for plain gzip and MANE preference) -> `results/r01_figure_20261102/
  promoterai/`: 3,798 rows in the window, 0 ref mismatches, 633/633 designed SNVs matched; six LDLR transcripts with
  identical scores (max |diff| 0); MANE ENST00000558518.6 (TSS chr19:11,089,462) used.
- Merged `evidence_v3`; `correlations_v2` (PromoterAI pairs prespecified in correlation_rules.json): vs MPRA 0.48
  (Sp1 site 0.85, outside 0.41); vs nucleosome-size -0.40; vs AlphaGenome RNA 0.06 (site 0.84). AUROC P/LP vs B/LB
  0.69. Notable: c.-99A>G and c.-98C>T (P/LP, CRE region) -0.69 / -0.77; c.-97G>C (LB) -0.77; SRE P/LP ~0.
- Figure A v9 (`results/r01_figure_20261102/fig_A_v9/`): panel b PromoterAI column filled: R1 65%, SRE 4%, R3 81%,
  F04 12%, CRE 67%, rest 10%.

## 2026-11-02 Data search for grant figures 1-2 (klone + Kopah)
- SLC39A4 (Fig 2): on klone, projects/SLC39A4_variant_slide/ (README; WT/MUT reps 1-3 BAMs, FiberHMM footprints,
  FIRE, analysis/site_stats.tsv: TF-size footprint over the variant WT 10.6/14.6/12.2% vs variant 2.7/3.1/2.4%;
  luciferase PNG from Fig. S15A, HEK293T, variant 21% of WT). Source reads: Kopah pEF_GFP_Gus archive.
- DPAGT1 / NAPA (Fig 1): not on klone; archived on Kopah (k_bmallo:bmallo-archive/ft_data/):
  plasmid_fiberseq/genomic_promoters (full-length promoters incl. DPAGT1, NAPA; 3 reps; merged per-construct BAMs
  kept, e.g. DPAGT1 2.1 MB), plasmid_fiberseq/IGVF_construct_mix (~300-bp constructs DPAGT11/13, NAPA1-3; 2 reps;
  source BAMs only, 0.9-4.3 GB), genomic/Hek293_Whole_Genome/target_bams (endogenous HEK293T reads at DPAGT1
  ["DPATG1"], NAPA, SLC39A4, etc.). Also plasmid_fiberseq/non_transfected_plasmid_controls (naked plasmids with and
  without Hia5, incl. SLC39A4WT/MUT; ~2 GB). Nothing restored yet.
- 2026-11-02 restore (user OK) to /gscratch/stergachislab/bmallo/large_home/data/kopah_restore/: genomic_promoters
  (DPAGT1, NAPA merged BAMs + reference; 15 files), Hek293_Whole_Genome (DPATG1/NAPA/SLC39A4 target BAMs; 13),
  IGVF_construct_mix (aligned source BAMs rep 1-2 + reference; 28), non_transfected_plasmid_controls (source BAMs
  both arms + reference + scripts; 11). rclone check --size-only: 0 differences for all four. Log
  data/kopah_restore/r01_restore_20261102.log.
- Read counts: genomic_promoters merged DPAGT1 488, NAPA 117; endogenous HEK293T full-coverage DPAGT1 19, NAPA 31,
  SLC39A4 23; IGVF ~300-bp constructs rep1/rep2: DPAGT11 31,675/67,634, DPAGT13 33,023/72,845, NAPA1 32,183/88,146,
  NAPA3 30,843/85,318, NAPA2 1,924/2,774. Naked controls: only contig pEFGFP (304,738 Hia5-treated, 546,943
  untreated); NO SLC39A4WT/MUT reads despite the README listing them as constructs (correction to the earlier plan).

## 2026-11-02 Figure 1 (DPAGT1) data
- Rules `tools/r01_figure/fig1_rules.json`. Inserts (`fig1_map.py` -> fig1/inserts.tsv, exact 25-mer, identity 1.0):
  full-length DPAGT1 chr11:119,100,225-119,104,342 (construct reversed); DPAGT11 (FIRE peak 1, gene body)
  119,100,271-570; DPAGT13 (FIRE peak 3, promoter) 119,102,160-459. DPAGT1 TSS (NM_001382.4, minus) 119,101,853.
- Motifs (`fig1_motifs.py`, P < 1e-5): CTCF 119,100,436-455 (peak 1, P 3.8e-6) and 119,102,636-666 (P 3.1e-7);
  GC box SP/KLF/MAZ 119,102,288-317 (peak 3).
- `fig1_fiberhmm.sh` (samtools -s seed 42, ~6,000 reads per construct per rep; fiberhmm-call hia5/pacbio, min-mapq 0,
  linear). `fig1_profiles.py`: kept DPAGT11 9,498/11,895 (not spanning 932, under-methylated 1,459, unresolved 6),
  DPAGT13 9,739/12,003, full-length 488/488, endogenous 29/30.
- Two bugs found by checks and fixed before any result was used: (1) MA spans are in original read orientation:
  reverse reads had 50-61% of m6A inside nucleosome calls vs 5-8% forward; after flipping 5-8% on both strands.
  (2) backbone positions of the 300-bp constructs were mapped onto genome coordinates; now masked to the insert.
  First version also ran for hours (per-base pandas lookups); vectorized (47 s).
- Descriptive (51-bp smoothed m6A, Spearman over positions, not independent): peak 1 endogenous vs 300-bp
  construct 0.98 (vs full-length construct 0.10); peak 3: 300-bp 0.58, full-length 0.57; full 4.1 kb vs full-length
  0.71. Endogenous depth only 22-25 molecules over each peak.
- Figure 1 draft `tools/r01_figure/fig1.py` -> `results/r01_figure_20261102/fig1_v2/` (v1 = layout draft). Smoothing
  changed to ratio of 51-bp sums (no gaps at G/C runs); peak 1 endogenous vs 300-bp rho 0.97 (0.98 with the
  earlier mean-of-fractions smoothing). Panel f: 44,458 footprints >= 10 bp (11,872 shorter pieces excluded).
  LDLR scale panel: 633/633 SNVs, median 1,955 molecules, 775 clones per SNV. No clear CTCF footprint at the peak-1
  CTCF motif in single molecules; panel a's TF/PIC labels are size classes, not identities.
- 2026-11-02 Figure 1 redesign (user): panel a = HPRC main-figure schematic row, cropped from
  HPRC_driver_variants/figure/main/panels_ab_paste_in_place.svg (`fig1_schematic.py`; labels "5,204 alleles" ->
  "WT + variants", "HEK293T" -> "cells", heading and the HPRC analysis arrow removed; vector schematic.svg kept for
  Illustrator) + a drawn final step; panel b = DPAGT1 gene model + proportion of molecules bound (nucleosome / PIC /
  TF size classes) for the endogenous genome (29 molecules) and the two 300-bp constructs (9,498 / 9,739), region
  chr11:119,100,225-119,102,550 (gene body in the data window + promoter peak 3). Old panels c-f dropped; script
  kept as `fig1_v1_draft.py`. Output `results/r01_figure_20261102/fig1_v3/`.
- 2026-11-02 Figure 1 draft 3 (user): one Illustrator-editable SVG (`fig1.py` -> fig1_v5/fig1.svg + PDF/PNG via
  rsvg). Schematic (`fig1_schematic.py` -> fig1/schematic2): HPRC vector icons; variant plasmids replaced by 7
  drawn reporter plasmids (coloured regulatory element + reporter-gene arrow); labels "Reporter / plasmid library",
  "Transfect / cells". Read stack: 14 DPAGT11 molecules (seed 20261102) from 2,139 eligible of 9,498 (span insert,
  nucleosomes over insert within +/-1 of median 2, >= 1 footprint 20-90 bp in the insert); FiberHMM footprints by
  LRSS bins, m6A, no FIRE. Panel b: LRSS profiles_DPAGT1 (paper style, 08_plot_tracks.py colours/smoothing),
  endogenous (n 23) and 300-bp peak 1 (n 74,010, all reads), window chr11:119,100,224-119,100,900.
- 2026-11-02 Figure 1 draft 4 (user): shared window = DPAGT1 gene body in the construct reads, chr11:119,100,224-
  119,101,853 (insert start .. NM_001382.4 TSS). Panel a reads from the 4.1-kb construct (14 of 163 eligible of 488;
  median 8 nucleosomes over the window), footprints only (no m6A, no FIRE). Panel b: footprint bins drawn x5 (labelled);
  genome track now also shows footprint bins (n = 24, so single calls are visible bumps). Output fig1_v7.

## 2026-11-02 — Figure 1 draft 5 (fig1_v8): clustered 4.1-kb reads, revised occupancy tracks
- Request: (1) cluster the 4.1-kb plasmid reads on the 300-bp region; (2) gene model above the reads with the 300-bp region highlighted; (3) footprint ×5 on the 300-bp track only; (4) solid fills, no lines, less smoothing, shortest peaks in front; (5) no n values, track labels styled like the paper figure.
- Rules were written to `fig1_rules.json` → `panel_a_clusters` before clustering was run. This is a display choice and is descriptive only.
- Inputs (md5: molecules.parquet 7831c060ecd35ecebf74d796bb41914a; profiles_DPAGT1.tsv.gz 0f3c4d9cad1774f1bd7203532869535c):
  - `fig1/profiles/DPAGT1_full.molecules.parquet`: 488 molecules.
  - `lrss2026_slides56/tables/profiles_DPAGT1.tsv.gz`: the `*_share` columns, re-smoothed here with Gaussian σ = 3 bp inside covered runs (LRSS used 10 bp).
- Clustering:
  - Region: 30 ten-bp bins over chr11:119,100,271–570; majority state per bin; one-hot encoded; Ward linkage; k = 4.
  - Clusters, in display order (molecules / share / shown / accessible-bin share; from `fig1_v8/clusters.tsv`):

    | Cluster | Molecules | Share | Shown | Accessible bins |
    |---|---|---|---|---|
    | 1 | 43 | 8.8% | 6 | 0.77 |
    | 2 | 165 | 33.8% | 20 | 0.33 |
    | 3 | 103 | 21.1% | 13 | 0.30 |
    | 4 | 177 | 36.3% | 22 | 0.09 |

  - 61 molecules shown in total. Each cluster gets round(60 × share), with a minimum of 6, so the total is 61, not 60.
  - Observation: sub-nucleosome footprints (20–90 bp) are the majority state in at most 6% of bins in any cluster. The clusters therefore separate accessible vs nucleosome-occupied 300-bp regions and nucleosome position. They do not separate TF-footprint states.
  - Cluster 1 (accessible) is where the 50–90 bp footprints inside the region sit.
- Tracks:
  - Layers are drawn in order of their maximum share, tallest first.
  - Footprint bins are ×5 on the 300bp track only; the legend note says so.
  - Labels: "Fiber-seq (genome-wide)"; "300bp region" with a "Plasmid Fiber-seq" bracket; y ticks 0/1 on the right; "Proportion of Reads Bound".
- Script: `tools/r01_figure/fig1.py` (git 78a4f18 + uncommitted). Output: `results/r01_figure_20261102/fig1_v8/` (svg, pdf, png, clusters.tsv, summary.json).
- An accidental draft-4 render into fig1_v8 was deleted before this run. It came from a Write that failed on a stale read.

## 2026-11-02 — Figure 1: footprint-size scan in the 300-bp region (exploratory, sets data-driven bins)
- Request: show only accessible reads with sub-nucleosome footprints, add the CTCF motif, and bin footprints by the sizes actually present.
- Script: `tools/r01_figure/fig1_fpsizes.py`. Output: `results/r01_figure_20261102/fig1/fp_sizes/`.
- Definitions:
  - Region: chr11:119,100,271–570.
  - A molecule must span the region (within 20 bp of each end).
  - "Accessible": nucleosome calls (≥90 bp) cover <50% of the region.
  - Sub-nucleosome footprint: <90 bp, with its midpoint in the region.
- Molecules (`summary.tsv`):

  | Dataset | Spanning the region | Accessible | Accessible with any sub-nucleosome footprint | Accessible with a footprint ≥16 bp |
  |---|---|---|---|---|
  | DPAGT1_full (4.1-kb construct) | 488 | 76 | 53 | 21 |
  | DPAGT11 (300-bp construct) | 9,498 | 5,787 | 4,586 | 2,622 |
  | Endogenous | 22 | 20 | 16 | 7 |

- Lengths (`length_hist_accessible.tsv`):
  - Most sub-nucleosome calls are under 16 bp. Median: 9 bp (4.1-kb construct), 10 bp (300-bp construct), 7.5 bp (endogenous).
  - Many recur at fixed positions with fixed lengths (`sites.tsv`). Likely gaps in AT-poor sequence rather than protein (not verified).
  - DPAGT11 modes ≥16 bp:
    - 19 bp (771 calls): mostly at midpoint 119,100,461–470, just 3′ of the CTCF motif.
    - Small peaks at 23 and 29 bp.
    - A low plateau over 31–47 bp.
    - A 48–57 bp mode (53 bp: 175 calls).
    - Sparse above 60 bp.
- CTCF motif (HOCOMOCO v14 CTCF.H14CORE.0.P.B, 119,100,436–455 −, P = 3.8e-6):
  - DPAGT11: 995 of 5,787 accessible molecules (17%) have a footprint ≥10 bp overlapping it by ≥10 bp. Median length 49 bp; IQR 24–54.
  - 4.1-kb construct: 2 of 76.
  - Endogenous: 2 of 20.
- Proposed data-derived bins (post hoc, pending user approval):
  - <16 bp: not shown (short, fixed-position calls).
  - 16–30 bp.
  - 31–60 bp (contains the CTCF ~50-bp mode).
  - 61–89 bp.
  - ≥90 bp: nucleosome.
- Implication: the 4.1-kb reads cannot show a CTCF footprint (2 molecules). The 300-bp construct reads can.

## 2026-11-02 — Figure 1 draft 6 (fig1_v9): user-chosen footprint bins, accessible 300-bp reads, CTCF motif
- Bins chosen by the user after the size scan (post hoc, data-derived): 15–30, 31–60, 61–89 bp, and ≥90 bp nucleosome. Calls <15 bp are not drawn or counted. Recorded in `fig1_rules.json` → `fp_bins_v2` and `panel_a_accessible` before rendering.
- Panel a reads: I switched to the 300-bp construct (DPAGT11), as recommended. The user set the bins without choosing a read source.
  - 2,783 of 9,498 molecules qualify: they span the region, are accessible (nucleosome calls cover <50% of it), and have ≥1 footprint of 15–89 bp with its midpoint in the region.
  - 40 are shown (seed 20261102), in Ward leaf order. Read IDs are in `panel_a_reads.tsv`.
  - The CTCF motif (HOCOMOCO v14, 119,100,436–455) is drawn as a track.
- Panel b tracks:
  - Footprint shares are recomputed from this project's FiberHMM calls (endogenous: 29 molecules, ≤24 per base; DPAGT11: 9,498).
  - FIRE comes from the LRSS table.
  - σ = 3 bp; sub-nucleosome bins ×5 on the 300bp track only.
  - Source data: `tracks.tsv`.
- 300-bp construct, raw shares (`tracks.tsv`):
  - 15–30 bp: maximum 0.166 at 119,100,468, just 3′ of the motif. This is the 19-bp mode; most of the purple calls in panel a sit here.
  - 31–60 bp: maximum 0.086 at 119,100,461; 0.081 at the motif centre (119,100,445).
  - FIRE: 0.52 over the motif.
- Trace: the 31–60 bp share at 119,100,445, recomputed directly from `DPAGT11.molecules.parquet`, is 765/9,498 = 0.0805. It matches `tracks.tsv`.
- Open question: is the 19-bp footprint at 461–470 part of the CTCF footprint (CTCF covers its motif plus flanks) or a separate factor? Not tested.

## 2026-11-02 — Figure 1 draft 7 (fig1_v10): CTCF-footprinted 300-bp construct reads along the whole plasmid
- Requests:
  - Panel a: CTCF-sized footprints only, 10–20 reads clustered on that footprint.
  - Endogenous track without sub-nucleosome footprints.
  - CTCF motif shown in panel b.
  - "300bp / region" stacked, with the bracket moved in.
  - Y label centred.
- The user first asked for 4.1-kb reads. `tools/r01_figure/fig1_ctcf_reads.py` (output `fig1/ctcf_reads/summary.tsv`) counts molecules with a 31–60 bp call overlapping the CTCF motif (119,100,436–455) by ≥10 bp:

  | Dataset | Molecules with a CTCF-sized footprint | Molecules spanning the motif | Share |
  |---|---|---|---|
  | 4.1-kb construct | 2 | 488 | 0.4% |
  | 300-bp construct (DPAGT11) | 765 | 9,498 | 8.1% |
  | Endogenous | 1 | 22 | — |

  The user then chose the 300-bp reads drawn along the whole plasmid. Rules: `fig1_rules.json` → `panel_a_ctcf` (revised before rendering).
- Panel a:
  - 15 of 765 molecules (seed 20261102), in Ward leaf order on per-base 31–60 bp occupancy over the motif ±60 bp.
  - Nucleosome and 31–60 bp calls only, drawn along the whole 4,774-bp linearized plasmid. Plasmid coordinate = hg38 − 119,098,913; insert 1,358–1,657; motif 1,523–1,542.
  - The footprint calls in `DPAGT11.molecules.parquet` already cover the backbone (only the profiles were limited to the insert).
  - The plasmid map shows backbone, insert and motif only, because the reference has no feature annotation. "NotI" at the ends is inferred from the construct name.
  - Read IDs: `fig1_v10/panel_a_reads.tsv` (not committed).
- Panel b:
  - Endogenous track: nucleosome and FIRE layers only.
  - CTCF motif drawn on the gene model and as a band across both tracks.
- Renders: PDF/PNG now use `tools/r01_figure/fonts.conf` (Arial → Arimo, metric-compatible). Earlier renders used DejaVu Sans, which is wider than Arial. That caused the off-centre y label and the legend crowding. The SVG itself was unaffected.

## 2026-11-02 — Figure 1 on hold (user: "not great"; pinned at draft 7, fig1_v10)

## 2026-11-02 — Figure 2 draft 2 (fig2_v2): SLC39A4 c.-169A>G, NF-Y footprint loss and luciferase
- Rules written before computing: `tools/r01_figure/fig2_rules.json`.
- Source: `projects/SLC39A4_variant_slide`.
  - HEK293T plasmid Fiber-seq, WT and variant co-transfected, 3 replicates.
  - Footprints called with FiberHMM `--circular`.
  - `occupancy.py` outputs are reused unchanged.
- Orientation QC (`fig2_orientation_qc.py`, rep 1, first 1,500 reads per BAM; `fig2/qc/orientation_qc.tsv`):
  - In these BAMs the MA spans are in stored-SEQ orientation. m6A inside nucleosome calls: 2.4–2.7% on both strands as stored; 25–27% on reverse reads if flipped.
  - So `occupancy.py`, which does not flip, is correct. The Figure 1 linear-mode BAMs needed the flip, so check per run.
- TF-sized footprint = 20–60 bp in every panel (the SLC39A4 README v5 definition). The slide's `site_stats.tsv` used 10–60 bp; this figure recomputes everything with 20–60.
- Panel d (`fig2_v2/site_stats.tsv`), reads with a 20–60 bp footprint over c.-169:

  | Replicate | WT | Variant | Odds ratio | Fisher p |
  |---|---|---|---|---|
  | Rep 1 | 162/2,130 (7.6%) | 36/2,241 (1.6%) | 5.0 | 8.6e-23 |
  | Rep 2 | 3,675/32,851 (11.2%) | 722/39,904 (1.8%) | 6.8 | <1e-300 |
  | Rep 3 | 16,900/181,155 (9.3%) | 3,039/213,396 (1.4%) | 7.1 | <1e-300 |

  - The test is within each replicate, with reads as observations. The 3 biological replicates are shown separately.
- Trace: rep 1 counts recomputed from `occ_rep_1_*.npy` at plasmid 6120 match exactly (162/2,130; 36/2,241).
- Panel c (`summary.json`; reps pooled, σ = 3 bp):
  - WT peak 9.5% at the NF-Y motif.
  - At the variant: WT 9.5% vs variant 1.6%.
- Panel e: WT − variant at the variant base is 5.8 / 9.1 / 7.8 percentage points (reps 1–3). The largest |Δ| outside the promoter is 2.2 / 0.9 / 0.7 percentage points; rep 1, with 2,130 reads, is the noisiest.
- Panel b, luciferase (Fig. S15A, HEK293T):
  - Supplemental.pdf page 27 is vector (378 drawings, no images). Exported with PyMuPDF 1.28.6 (installed with `pip --target` into the session scratchpad, not into any conda env) by `fig2_s15_svg.py`.
  - Each per-glyph text run is split into word-level tspans so it renders correctly in rsvg and in Illustrator.
  - Embedded clipped to x 56–196 / y 68–284 pt, left of the HEK293T/HepG2 divider. One closing line was added on the right frame edge.
  - The HepG2 half (p = 0.056) is not shown.
  - The raw luciferase values are not on klone, so "21% of WT" from the slide is not reproduced.
- Input md5: reads.tsv.gz and Supplemental.pdf are listed on the line below.
  - 9af2f0d27f340b8771795d3481178986 /gscratch/stergachislab/bmallo/large_home/projects/SLC39A4_variant_slide/analysis/reads.tsv.gz
  - bb501840e8f9a6b46d4ef9dbafb1aab3 /gscratch/stergachislab/bmallo/large_home/projects/SLC39A4_variant_slide/source/Supplemental.pdf

## 2026-11-02 — Figure 2 draft 3 (fig2_v3)
- Footprint-size question (user item 6), from `tools/r01_figure/fig2_nfy_size.py` (`fig2/nfy_size/`):
  - The WT-specific footprint over c.-169 has no single size. The WT excess is spiky and spread over 3–58 bp, because FiberHMM calls the protected stretch as fragments of variable length.
  - The excess stops at ~60 bp. Above that (61–120 bp) it is +1.0 percentage points (pp) at 1.05×.
  - Pooled excess by bin: 1–60 bp 14.4 pp (5.3×); 20–60 bp 8.1 pp (6.5×); 25–50 bp 5.9 pp (7.0×); 1–9 bp 4.4 pp (6.1×).
  - Decision: keep the pre-specified 20–60 bp. A narrower bin barely raises the fold and drops real signal.
- Panel a: SnapGene-style map, to scale (plasmid 6,400→2,800 in gene orientation) with thin outlined feature arrows. The motif is labelled only "NF-Y motif".
- Panel b:
  - The side text is removed.
  - The closing frame edge sits at the S15A frame rectangle (y 80.07–252.83 PDF pt, 1 pt stroke). In draft 2 it was offset.
- Panel c: read counts removed from the legend.
- Panel e: the mean of the 3 per-replicate tracks (replicates weighted equally). Values are on the line below.
  - Mean Δ at plasmid 6120: 7.55 percentage points. Largest |mean Δ| outside the promoter: 1.14 percentage points.

## 2026-11-02 — Figure 2 draft 4 (fig2_v4)
- Removed the draft stamp and the note under panel e (the user asked). Page height trimmed to 5.5 in. Commit and date stay in summary.json; no numbers changed.

## 2026-11-02 — Figure 3 draft 1 (fig3_v1): LDLR promoter truth-set validation
- Rules: `tools/r01_figure/fig3_rules.json`, written before running. Script: `tools/r01_figure/fig3.py`.
- No new statistics. The figure redraws the reviewed LDLR panels from their saved source tables in the Figure 2 style (editable SVG, Arial; renders via fonts.conf).
- Panels:
  - a: HepG2 MPRA along the promoter, ClinVar class by marker. 633 SNVs: 529 not in ClinVar, 84 VUS, 6 conflicting, 6 B/LB, 8 P/LP.
  - b: region × readout share in each readout's top 20%. Post hoc / exploratory per `figA_rules.json`.
  - c: Sp1 site, 27 SNVs, Fiber-seq protection vs MPRA. Spearman ρ = −0.69, P = 7.0e-05 (`correlations_v2`).
  - d: mechanism category by ClinVar class (`panel_b_v2`, |log2| ≥ 0.5). Repressor categories are labelled as candidates.
  - e: HepG2 vs HEK293T local inaccessibility, 113 SNVs with ≥3 HEK clones. ρ = 0.62 (95% CI 0.47–0.75). Separate axes, not pooled.
- c.-156C>T is not used as a repressor example (likely Hia5 sequence artifact; see earlier entries).
- Input md5s are on the line below.
  - fd3fe0258a8a5f5f687f0d97553dc49f fig_A_v9/panel_a.tsv; 999974d5a28fc696fe8f385afa092d8a fig_A_v9/panel_b_regions.tsv; 87c2f22c72f8d05ef8173eba68d7380f fig_A_v9/panel_c_site.tsv; 1b2e2b40120d5095339c2d503713f488 fig_A_v9/panel_d_hek.tsv; 71e2454de68c63b4a24009e53d043c91 panel_b_v2/panel_b_counts.tsv; cf8395ea6b6655d28450f7c6d904e56f correlations_v2/pairs.tsv; 

## 2026-11-02 — Figure 3 draft 2 (fig3_v2): rule-based regions, matched in panels a and b
- User request: every region measured in b must be highlighted in a, easy to define, and contain significant ClinVar variants.
- Rule (`fig3_rules.json` → `regions_draft2`, written before computing): a WT motif-scan match at P < 1e-5 of an LDLR-promoter TF family, containing ≥1 ClinVar P/LP SNV.
- Kept:
  - SRE-1 (SREBP motif), oligo 127–136: P/LP 133, 134, 136.
  - GC box (SP1, Repeat 3), oligo 143–151: P/LP 149, 150; conflicting 145, 146.
- Dropped:
  - SP-like 89–100: no P/LP.
  - "CRE" 182–192: P/LP 186/187, but no motif at P < 1e-5.
  - F04 core: no motif, no P/LP.
- Panel b recomputed in `fig3.py` with fig_A.py's readout definitions; output `fig3_v2/panel_b_regions.tsv`. Share of each region's SNVs in each readout's top 20%:

  | Region | SNVs | MPRA | Fiber-seq | AlphaGenome | PromoterAI |
  |---|---|---|---|---|---|
  | SRE-1 | 30 | 77% | 13% | 20% | 3% |
  | GC box | 27 | 78% | 81% | 78% | 81% |
  | Rest of element | 576 | 15% | 18% | 17% | 18% |

- Check: the GC box row equals fig_A_v9's R3 row exactly, so the recomputation matches.

## 2026-11-02 — Figure 3 draft 3 (fig3_v3): all three motif sites
- User: keep all three motifs. Panels a and b both use the SP-like site (89–100, ZFP3 match P 3.4e-6), SRE-1 (127–136) and the GC box (143–151). Rules: `fig3_rules.json` → `regions_draft3`.
- Panel b (`fig3_v3/panel_b_regions.tsv`), share of each region's SNVs in each readout's top 20%:

  | Region | SNVs | MPRA | Fiber-seq | AlphaGenome | PromoterAI |
  |---|---|---|---|---|---|
  | SP-like | 36 | 69% | 28% | 64% | 78% |
  | SRE-1 | 30 | 77% | 13% | 20% | 3% |
  | GC box | 27 | 78% | 81% | 78% | 81% |
  | Rest of element | 540 | 11% | 17% | 14% | 14% |

- F04 and the CRE are not drawn.
  - The CRE is the analysis_v2 "CREB1/CREM 182–192" window (c.-103..-93), sequence AAATGACGTGGG. Its TGACG core holds ClinVar P/LP c.-99A>G and c.-98C>T. In this scan the best matches are XBP1 (P 2.7e-5) and CREB1 (P 7e-5).
  - F04 is oligo 66–82 (c.-219..-203), a data-defined footprint with no motif and no P/LP.

## 2026-11-02 — Figure 3 draft 4 (fig3_v4): F04 and CRE added
- User: "Add both". Five regions now, 5′→3′, shaded in a and measured in b. Rules: `fig3_rules.json` → `regions_draft4`.
- F04 is drawn with a dashed box because it is data-defined. The TSS label moved to the right of the arrow to clear the CRE label.
- Panel b (`fig3_v4/panel_b_regions.tsv`), share of each region's SNVs in each readout's top 20%:

  | Region | SNVs | MPRA | Fiber-seq | AlphaGenome | PromoterAI |
  |---|---|---|---|---|---|
  | F04 footprint | 51 | 24% | 67% | 10% | 12% |
  | SP-like | 36 | 69% | 28% | 64% | 78% |
  | SRE-1 | 30 | 77% | 13% | 20% | 3% |
  | GC box | 27 | 78% | 81% | 78% | 81% |
  | CRE | 33 | 0% | 6% | 64% | 67% |
  | Rest of element | 456 | 10% | 12% | 11% | 11% |

- F04 and CRE match fig_A_v9 exactly.

## 2026-11-02 — Figure 3 draft 5 (fig3_v5): truth-set detection matrix and element-wide MPRA × Fiber-seq map
- Why: the user asked to rethink panels a/b against the grant. The grant's multi-modal claims are about known variants:
  - Barrier 2: single readouts miss variants.
  - Aim 1B: proxy benign classes.
  - Aim 1C: per-assay false-negative rate, discordance, joint calibration.
  - Aim 3A: comparison with the predictors.
  The approved plan is at `~/.claude/plans/calm-popping-canyon.md`. User decisions: matched top-20% flag; P/LP and B/LB columns only; drop the mechanism bars; keep HEK293T small.
- Rules: `fig3_rules.json` → `layout_draft5`, written before computing. Post hoc display choice: the region view and a read-only preview of control counts had been seen.
- Layout:
  - a: unchanged.
  - b: control-variant detection matrix.
  - c: element-wide MPRA × Fiber-seq map.
  - d: HepG2 vs HEK293T (unchanged, formerly e).
  - The mechanism bars and the GC-box-only scatter leave the main figure.
- Flag rule: |score| ≥ the readout's 0.8 quantile across 633 SNVs. Cutoffs: MPRA 0.97 log2; Fiber-seq 2.44 z; AlphaGenome 1.22 SD; PromoterAI 0.228.
- Panel b (`fig3_v5/panel_b_matrix.tsv`, `summary.json`):

  | Readout | P/LP flagged | B/LB flagged | gnomAD-observed flagged | All SNVs flagged |
  |---|---|---|---|---|
  | MPRA | 4/8 | 1/6 | 23% | 20% |
  | Plasmid Fiber-seq | 4/8 | 2/6 | 25% | 20% |
  | AlphaGenome | 4/8 | 0/6 | 21% | 20% |
  | PromoterAI | 4/8 | 2/6 | 28% | 20% |
  | Any readout | 7/8 | 3/6 | 50% | 46% |

  - The gnomAD-observed set is 107 SNVs, all rare.
  - c.-151C>G (SRE-1, MPRA −0.48) is flagged by no readout.
  - Benign flags: c.-217C>T (MPRA ▲, Fiber-seq ▼, PromoterAI ▲); c.-205G>T (Fiber-seq ▼); c.-97G>C (PromoterAI ▼).
- Panel c, flagged-by counts: both 46, MPRA only 82, Fiber-seq only 81, neither 424.
- Verification: an independent recompute from `evidence_v3` + fp_features reproduces every count, with 0/70 matrix cells different. Trace: c.-136C>T has MPRA |−2.62| ≥ 0.97, so it is flagged.
