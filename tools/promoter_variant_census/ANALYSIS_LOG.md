# Promoter variant census: analysis log

Question (user, R01 planning): for one 300 bp site-saturation window per gene (LDLR, HBB, GATA2, SLC39A4), how many
variants are in gnomAD, at AF > 1%, not GTEx eQTLs, and predicted "no effect" by AlphaGenome and PromoterAI.

## 2026-11-04 Step 1: windows
- TSS/strand from Ensembl REST release 116 (GRCh38) MANE Select; first rule −250..+49 around the TSS (TSS = 0), no
  coding bases. Kircher 2019 element coordinates from the satMutMPRA portal (LDLR construct chr19:11,089,231–548;
  HBB chr11:5,227,022–208).
- GATA2: IS = NM_032638.5 (MANE Select, distal, 128,493,201); IG = NM_001145661.2 (MANE Plus Clinical, 128,488,530).
- User supplied FIRE peaks; windows shifted to cover them (FIRE coordinates taken as 1-based inclusive):
  LDLR → 11,089,249–548 (FIRE 287/297 bp; last 10 FIRE bp are coding); GATA2 IS centred on FIRE (300/316);
  GATA2 IG centred (211/211); SLC39A4 two windows on upstream peak B (W1 −250..+49, W2 −492..−193, 58 bp overlap),
  coding peak A not included; HBB unchanged (no FIRE peak). Locked in windows.tsv after user OK.

## 2026-11-04 Step 2: definitions
- rules.json written before any variant in the windows was retrieved; approved by the user the same day.
- Checks made before writing: gnomAD GraphQL fields validated on an unrelated region (chr1:1,000,000); GTEx API v2
  variant table carries `maf01`; V10 fine-mapping endpoint returns nothing (V8 only); AlphaGenome biosample metadata
  listed (no adult erythroid sample; HEK293 only in CAGE); PromoterAI tss.tsv contains all six window transcripts.

## 2026-11-04 Step 3: counts (implementation notes)
- fetch.py gnomad: six region queries padded ±50 bp, raw JSON kept in results/.../gnomad/.
- fetch.py gtex: significant pairs per target gene (V10, V8) and V8 fine-mapping per gene in one paged query each;
  variant-table lookups (V10, V8) only for records with genome AC ≥ 1 (the census uses genomes only; exome-only
  variants are context and need no GTEx status); per-variant any-gene eQTLs only when present in the V10 table.
- fetch.py promoterai: awk filter of promoterAI_tss500.tsv.gz to the windows (14,730 rows; 900 SNVs for every window
  transcript, 1,626 for SLC39A4 W1∪W2).
- ag_score.py: all possible SNVs + gnomAD indels with genome AC ≥ 1. SLC39A4 W1 and W2 share one 1 Mb interval
  (centre of W1∪W2) so overlapping SNVs are scored once; other windows each have their own interval centred on the
  window. Saves kept-biosample tracks and per-output-type max |quantile| over all tracks.
- GTEx variant_lookup.tsv header bug (header taken from a row without the eQTL columns): fetch.py fixed; the file was
  rebuilt from variant_lookup.raw.tsv without re-querying (644 rows). Gene-level and per-variant significance agree
  for all 517 gnomAD-PASS variants (2 significant).
- AlphaGenome: 5,259 variant-contexts, 5,867 s, 0 failed requests.
- assemble.py -> results/promoter_variant_census_20261104/ (counts.tsv md5 43962f46..., funnel.tsv da8646a9...,
  variants.tsv d8aa8c67..., all_snvs.tsv 9ef7e9c5...). Checks: 0 ref mismatches (gnomAD, PromoterAI), 0 missing scores.

## 2026-11-04 Results (prespecified definitions)
- gnomAD v4.1 genome PASS SNVs / 900: LDLR 61, HBB 77, GATA2 IS 107, GATA2 IG 76, SLC39A4 W1∪W2 141/1626.
  Indels 2, 2, 11, 11, 6. AF > 1% overall: 0, 0, 0, 1, 2 SNVs; any group (AN >= 2000): 2, 2, 2, 2, 6. No common indel.
- GTEx V10: 2 SNVs significant for the target gene (GATA2 IG rs73862209; SLC39A4 rs117568226), 0 tested-not-significant,
  all other gnomAD variants not tested (only 3 in-region variants carry maf01).
- AlphaGenome no effect (|q| < 0.99 in every gene-set RNA/CAGE track), all SNVs: LDLR 33, HBB 11, GATA2 IS 382,
  IG 194, SLC39A4 543/1626. Driven by RNA_SEQ for LDLR/HBB/SLC39A4 W1 (median max |q| 0.999-1.000).
- PromoterAI |score| < 0.1, all SNVs: 548, 762, 576, 739, 1454/1626.
- Intersection (common, not eQTL, both no effect): 0 in every window with overall AF; 3 with any-group AF
  (GATA2 IS rs563914744; SLC39A4 W2 rs142890446, rs145429748); 0 with the strict (tested) eQTL rule.

## 2026-11-04 Truth-set potential (truth_set_rules.json, approved "Approve all")
- User framing: candidate benign = common (gnomAD) + not a GTEx eQTL + no predicted RNA/DNase change (AlphaGenome,
  PromoterAI); candidate pathogenic = ClinVar P/LP. B3 replaces the census AlphaGenome rule post hoc (quantile 0.99
  passed almost no promoter SNV); top-20%-within-element rule as Figure 3 panel b; DNase replaces CAGE.
- fetch.py clinvar (E-utilities, 2026-11-04): LDLR 161 records, HBB 102, GATA2 IS 4, IG 0, SLC39A4 W1 5, W2 0.
- truth_sets.py -> results/promoter_variant_census_20261104/truth_sets/ (truth_sets.tsv md5 d5b963dc...). Checks: 0 SNVs missing AlphaGenome or PromoterAI.
- P/LP (SNV+indel): LDLR 14 (10+4), HBB 29, GATA2 IS 0, IG 0, SLC39A4 0 (Figure 2 c.-169A>G is not in ClinVar).
  B/LB: 8, 5, 1, 0, 2. VUS+conflicting: 130, 59, 2, 0, 3.
- Candidate benign SNVs: common (any group) 2, 2, 2, 2, 6 -> not eQTL (any gene) 2, 2, 2, 1, 5 -> predicted no effect
  1, 1, 2, 0, 3 (7 total; 3 already ClinVar B/LB). Strict (GTEx-tested) and overall-AF versions: 0 everywhere.
  12/14 common SNVs untestable in GTEx; the 2 tested are eQTLs (GATA2 IG rs73862209; SLC39A4 rs117568226).
- LDLR c.-217C>T (rs17249141, ClinVar benign, expert panel) is the only measured control meeting B1; it fails B3
  (AlphaGenome DNase large, PromoterAI +0.333) and has an MPRA increase.

## 2026-11-04 GATA2 +9.5 intronic enhancer replaces the GATA2 promoters (truth_set_rules.json gata2_enhancer)
- Located from ClinVar P/LP c.1017+572C>T (chr3:128,483,288) and c.1017+513_540del (NM_032638.5; c.1017 = exon-4 end
  128,483,860); E-box CATCTG at c.1017+530..535, GATA AGATAA at +544..549; ENCODE dELS EH38E2236526. User FIRE peak
  chr3:128,482,922-128,483,424 -> GATA2_E95_W1 128,482,922-221 and W2 128,483,125-424 (97-bp overlap).
- fetch.py gained --windows (gnomad, clinvar) and reuses the gene-level GTEx tables. Bug: the 177 appended GTEx rows were
  written in dict order, not the header order (misaligned columns; earlier elements' b1_untested became 0 in a first
  rerun). File realigned from variant_lookup.raw2.tsv (no re-query); fetch.py now writes LOOKUP_COLS order. After the
  fix all earlier-element counts equal the previous run.
- AlphaGenome: 1,523 new variant-contexts, 1,558 s, 0 failures. PromoterAI not applicable (counted as not failing).
- Enhancer (W1+W2, 1,509 possible SNVs): ClinVar P/LP 2 (1 SNV + 1 indel), B/LB 85 (79 SNVs + 6 indels, all single
  submitter, GATA2 deficiency), VUS 5, conflicting 9. 31 B/LB SNVs lie in the E-box/GATA/ETS core (c.1017+512..+579),
  incl. +533, +534, +535 in the E-box; 35/79 B/LB SNVs have a top-20% AlphaGenome change (29% of all SNVs).
  gnomAD 128 SNVs, 14 indels; common 2 (rs62270945 overall 1.77%, GTEx eQTL for IFT122; rs116644718 afr 1.05%, not
  tested) -> candidate benign 1 (rs116644718). truth_sets.tsv md5 4c52dcc2...
