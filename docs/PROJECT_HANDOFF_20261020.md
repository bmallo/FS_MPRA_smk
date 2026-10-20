# LDLR plasmid Fiber-seq / MPRA: scientific and implementation handoff

**Status through 20 October 2026. Read this first.** This is the current synthesis; dated reports preserve earlier decisions and may propose steps that were subsequently completed. The code and reports distinguish observed differences from statistically supported discoveries. **No replicate-3 Fiber-seq variant effects have been released as statistically significant.** Published functional measurements retain their source statistics; those do not validate our footprint effects.

The most useful biological observations are a local protection/detection change at **chr19:11089389** and a broader accessibility-associated region at **chr19:11089407–11089414**. Several internal checks support these candidates, but sequence-dependent labeling/calling and single-preparation inference remain unresolved. The next priority is focused experimental validation, not filter tuning until significance appears.

## 1. Scientific question and experimental facts

The objective is to determine how individual base substitutions alter small, TF-sized protection and overall promoter accessibility, then relate architecture to functional activity and candidate binding motifs. Footprint coverage is an assay measurement, not instantaneous occupancy of an identified protein.

- Project data: `/gscratch/stergachislab/bmallo/large_home/projects/LDLR_MPRA`; highest-priority preparation: `fiberseq/HepG2/replicate_3`.
- HEK293T uses an older LDLR library; HepG2 uses a newer, shorter LDLR element. Do not combine their coordinates or barcode maps without explicit library reconciliation.
- HepG2 replicates 1–3 are independent biological preparations. Runs 1–2 used **10 minutes, 1.5 µL Hia5 per two million cells**; run 3 used **5 minutes, 3 µL per two million cells**. Enzyme stock activity and remaining preparation equivalence are incompletely documented. Replicate 3 is not a condition-matched third replicate for pooled inference.
- The primary interval is the user-confirmed **211-bp mutable element**, oligo positions 21–231. Multiple constructed barcodes per variant should be readable on long molecules.
- High-methylation molecules are from the **same Fiber-seq reaction**, not separately prepared naked-DNA controls. Their cellular compartment is unknown. Extranuclear/endosomal plasmids are a hypothesis, not an established classification.
- Transient binding during the five-minute reaction, uneven labeling efficiency, heterogeneous plasmid chromatinization, sequence-dependent A/T observability, sequencing/calling errors, read recovery, barcode representation and shared-cell dependence can all affect comparisons. Neither a statistical test nor filtering can recover erased footprints or identify these latent causes by itself.
- Only existing published functional data are available. `mpra_published/GRCh38_LDLR.tsv` and `GRCh38_LDLR.2.csv` are retained separately and are not assumed independent or matched biological samples. Future user-generated RNA/DNA data should use MPRAsnakeflow; that pipeline has not been run here.

## 2. Work completed, in order

| Phase | What changed / what was learned | Authoritative report |
|---|---|---|
| Repository review and pilot | Separated inherited exploratory scripts from a new sequence catalog, deterministic molecule sampling, explicit missingness, and local dashboard. Corrected rep3 1% pilot: 28,415 measured molecules. | [Scope](redesign_scope_and_plan.md), [pilot](pilot_run_20261014.md) |
| Early run comparison | Bounded SQLite imports processed full run-1/2 ledgers: 210,294 and 21,078 measured WT/SNV molecules. Run 2 was dominated by short fragments; full-plasmid recovery 2.8% versus 33.6% in run 1. Regional agreement was poor. User redirected effort to rep3. | [Full-data/recovery audit](full_data_20261016.md), [condition comparison](condition_comparison_20261014.md) |
| Raw-BAM processing foundation | Added separate Snakemake workflow for preflight, catalog, alignment, QC/genotype/barcode extraction, caller execution and feature export. Real 10,000-record smoke test completed. | [Processing workflow](processing_workflow.md) |
| Barcode benchmark | Compared exact grouping, quality-directional correction and pinned Pacybara clustering core. Kept exact grouping after yield/false-merge audits. | [Benchmark findings](barcode_benchmark_20261016.md) |
| Complete rep3 ledger and footprint join | Reused existing alignment and historical FiberHMM calls; preserved all primary molecules and original IDs. No full-data realignment or recall was performed. | [Ledger](rep3_ledger_20261016.md), [join](rep3_footprints_20261017.md) |
| Initial statistical engine | Built reusable barcode arrays, fixed endpoints, Welch intervals/tests, BY correction and release gates. Calibration failed; candidate significance withheld. | [First statistics](rep3_statistics_20261017.md) |
| Heterogeneity and caller audit | Characterized zero masses, methylation/protection distributions, barcode variation, call lengths, missingness and score provenance. | [Heterogeneity](rep3_heterogeneity_20261017.md), [quality](rep3_footprint_quality_20261018.md), [caller provenance](rep3_caller_provenance_20261018.md) |
| Expanded calibration | Compared Welch/studentized permutation, presence/coverage, realistic imbalance, sparse nulls, small effects and interval coverage. No general inference release. | [Expanded benchmark](rep3_heldout_calibration_20261018.md) |
| Biological landscape and atlas | Categorized WT architecture, all 633 SNVs, candidate regions and functional associations. Added previous/next navigation and synchronized functional lollipops. | [Landscape](rep3_landscape_20261019.md), [navigation](landscape_navigation_20261019.md) |
| Genomic/motif/barcode audits | Independently verified sequence orientation and genomic joins; added fixed-version motif scans and barcode/sampling/observability sensitivities. | [Coordinates](genomic_coordinate_audit_20261020.md), [biology](rep3_candidate_biology_20261020.md) |
| Raw-call candidate audit | Extracted 134,529 molecules for 27 predefined candidate SNVs plus WT samples; compared unchanged A/T sites and 18 cohort/weighting settings. | [Latest audit](rep3_raw_candidates_20261020.md) |

The inherited `workflow/Snakefile`, `workflow/scripts/`, older motif/co-occupancy tools and legacy statistical documents remain for history. Their earlier calibration claims were **not revalidated for the new barcode-aware cohort**. The new workflows are the explicitly named `.smk` files below. Do not treat inherited output as another independently reproduced result.

## 3. Inputs, identity, coordinates and cohort

The physical reference `LDLR_MPRA_NotI` is 4,647 bp. The mutable interval is plasmid **[1234,1445)**; the N15 barcode is **[3294,3309)**, reference-forward. The catalog contains 1 WT, 633 SNVs, 146 deletions and 637 insertions (1,417 unique sequences). Stable IDs encode library identity and a sequence SHA256; oligo aliases are retained separately. Primary biology currently covers WT/SNVs, not designed indels.

The mutable sequence exactly matches **GRCh38 chr19:11089285–11089495, forward strand**, verified against plasmid/catalog/local genome and an independently retrieved UCSC response including flanks. All 633 SNVs in each published source reconcile by reference allele and altered sequence.

```
genomic_pos1 = 11089264 + oligo_pos1
genomic_pos1 = 11089285 + element_relative_pos0
element BED = chr19 11089284 11089495
```

Displayed genomic coordinates are 1-based inclusive; array/BED intervals are 0-based half-open. The legacy `extras/data/LDLR_mpra_scores.tsv` instead requires `11092732 - legacy_position` and allele complementation; all 954 score/p-value pairs reconcile with source 2. **Complementing alleles does not negate activity.** Never assign genomic coordinates to plasmid backbone/barcode bases.

The complete ledger counted **8,715,827 primary plasmid alignments**, excluding 13,540 secondary/supplementary alignments. Genotype observations include 4,739,908 exact designs with adequate element base quality, 1,325,305 exact designs with low quality, 2,344,144 off-design observations, and 306,470 unresolved flanks. Off-design observations are not automatically genuine unintended constructs.

| Analysis population | Molecules | Supported barcodes |
|---|---:|---:|
| Sequence-qualified exact WT/SNV, broader heterogeneity cohort | 2,660,013 | Not required |
| Primary eligible WT/SNV, barcode support not required | 2,594,368 | Not required |
| Primary eligible supported-barcode cohort | 2,102,906 | 568,354 |
| WT within supported cohort | 985,264 | 268,283 |
| SNVs within supported cohort | 1,117,642 | 300,071 |

All 633 SNVs have eligible molecules; 632 have supported barcode data; 630 meet the 50-barcode/100-molecule display gate. **160A>G** has no supported eligible molecules, **222T>C** has one barcode/two molecules, and **176G>A** has four barcodes/eight molecules. Missing measurements are not zero effects.

The footprint join matched every primary read ID, with no sequence, barcode, length or strand disagreements. Of all primary rows, 3,707 lack MA annotations, all outside the primary eligible cohort. A present evaluated-empty track is zero; absent annotations are missing. Call lengths are classified **before clipping** to the element; circular splitting and clipped boundaries limit physical interpretation.

## 4. Filtering, barcode assignment and normalization

Primary sequencing requirements: exact designed WT/SNV genotype, mutable/barcode baseQ ≥20 where required for assignment, polymerase passes ≥10, read quality ≥0.99, mapQ ≥20, circular/full-plasmid alignment/length criteria, known methylation state and complete element for footprint analysis. Structural settings include 50-bp end tolerance, ≤100-bp soft clipping and concatemer threshold 1.1× reference length; consult the extraction code/config for exact Boolean logic. The configured methylation range is **0–1**, not an aggressive biological methylation cutoff. Accepted barcode assignment additionally requires at least two qualifying observations and no conflicting qualified genotype. Unresolved genotype is never rescued solely from barcode identity.

Barcode benchmarking used an unbiased 1% read sample (87,618 molecules), targeted full-depth neighborhoods around 400 discovery barcodes, held-out linkage, and 2,980 synthetic molecules. The random sample was mostly singleton barcodes. In the targeted cohort, Pacybara added 23 assignments but changed 18 barcode observations, including merges of already supported neighboring families. It did not improve held-out disagreement. Synthetic tests demonstrated real-family merging and ambiguous-parent errors. Quality-directional correction added no real-data yield in these cohorts. **Exact grouping remains the default**; only the pinned Pacybara core/settings were benchmarked, not the entire upstream workflow or every possible correction approach.

Full-ledger supported assignments total 4,529,384 molecules across all retained classes; the headerless MPRAsnakeflow export has 986,680 barcode→element mappings. This differs intentionally from the smaller eligible WT/SNV cohort. `catalog/design.fa` uses matching element IDs. Confirm future RNA/DNA barcode orientation and physical library before importing; barcode absence in Fiber-seq is not functional absence. Functional count import remains unimplemented.

**Normalization:** average eligible molecules within each barcode, then average barcodes equally. Coverage is covered bases / region bases; presence is a per-molecule indicator averaged within barcode. Alternate equal-molecule summaries are explicit. Fractions are normalized for usable support, not raw read counts. Barcodes and molecules are not independent biological preparations; shared-cell dependence is unmeasured. Small footprints/MSPs can overlap and do not partition the element into mutually exclusive states.

## 5. What the raw data and annotations reveal

- WT small-call (10–60 bp) mean element coverage is **2.60%**, with **21.84%** of weighted molecules having any element call. MSP coverage averages **52.59%**; 20.97% have none. Sparse local TF signal is the principal statistical difficulty.
- WT nucleosome-call counts have 10th/median/90th percentiles **0/8/17**. **14.57%** have no ≥85-bp nucleosome calls anywhere, yet **29.79% of that group** has an element small footprint. “No nucleosomes” is not a validated protein-free control.
- **5.93%** have exactly one nucleosome-class segment longer than 1,000 bp, with near-zero distal labeling and near-complete distal protection. This is unresolved protection, not one histone octamer. Across the element, nucleosome-sized 85–250-bp coverage is **19.75%**, while >1,000-bp coverage is **13.66%**.
- Heterogeneity persists within barcode families: 97.54% vary in nucleosome-call count, 99.94% in distal methylation, and 91.28% in MSP coverage. Median WT barcode depth is three. These spreads cannot separate technical from biological variation.
- Exact whole-plasmid coverage would be a destructive filter: only 31 supported eligible molecules cover every reference base, whereas 99.886% cover ≥99%. Regional measurements use actual aligned opportunities.
- A 10,000-WT-barcode sample gives A/B TF-profile Pearson **0.9987**, but both halves share sequence, preparation and caller. TF AQ≥100 lowers element coverage from **2.61% to 1.66%**; that threshold is not biologically calibrated. All 331,304 sampled nucleosome-class calls have AQ first byte zero, plausibly an unverified-score sentinel; MSP has no equivalent score.
- Target m6A ML bytes are restricted to **240–255** in bounded samples of available raw-ft, aligned and footprint BAMs, and in all positive element calls in the latest 134,529-molecule candidate extraction. ML≥128 versus ≥204 is therefore uninformative here. The restriction precedes alignment/calling in available files. Earlier probabilities/kinetics and exact upstream settings have not been recovered.

Never use `pysam.modified_bases` on these mixed tags in this environment: it can segfault. The custom parser consumes all MM groups and matching ML columns, handles reverse orientation, and preserves missing channels/unknown skipped sites. A known site without a retained call is not a measured original probability of zero.

## 6. Statistical analyses already attempted

| Analysis | Design and result | Consequence |
|---|---|---|
| Initial barcode Welch + BY | 21 fixed overlapping 20-bp TF coverage windows + whole-element MSP; 633×22 = **13,926** hypotheses. Coverage/sparse-support gates, pointwise intervals, BY including untestable members. Candidate q≤.05 in 330 rows/174 variants. WT A/A and hierarchical calibration failed: 40/88 WT strata failed and 19 were insufficiently testable; worst observed rejection 13% at nominal 5%. | **All significance flags remain zero. The 330 rows are not discoveries.** |
| WT quality benchmark | 7,200 contrasts, two barcode-disjoint halves of 10,000 WT barcodes; Welch versus studentized permutation (499 permutations), coverage/presence and heterogeneous simulations. Several inflated null rates improved; worst observed permutation WT rate 8.5%. | Candidate method, not general validation; p floor .002 inadequate for the large family. |
| Expanded WT calibration | 258,283 WT barcodes/948,189 molecules excluding the 10k sample; 4,000 splits ×43 outcomes =172,000 correlated contrasts; 50/100/400/1,000 pseudo-variant barcodes vs entire complement, 9,999-draw reference banks. Also 32 simulated nulls +18 alternatives, 200 draws each, 999 permutations. | More realistic imbalance and sparse regimes; disjoint from recent quality sample, **not globally untouched WT**. |
| Expanded error/power/interval checks | Welch WT rejection up to **26.4%** and sparse synthetic null up to **68.5%**; permutation greatly improved but some elevated/conservative regimes remained. MSP permutation rates 4.5/5.7/4.9/5.3%. A 2%→1.5% effect detected in 3.5/22.5/40% at 100/400/1,000 barcodes, uncorrected .05. | Sparse variant support cannot be fixed by adding WT reads. Nonsignificance is not equivalence. |
| Interval alternatives | Welch coverage as low as 73.6% in WT and 31.5% in synthetic nulls. Hoeffding intervals covered every tested truth but were too wide (roughly 10–42 pp across WT group sizes). | Better permutation p-values do not repair Welch intervals. No generally useful calibrated interval released. |

Studentizing permutation does not make unequal-distribution mean testing finite-sample exact. Shared finite reference-bank uncertainty, direction-specific power, extreme multiple-testing tails, and biological dependence remain unresolved. Fresh-bank checks, bootstrap/inverted-test intervals and calibrated hierarchical/hurdle models were proposed but **not established or run as completed alternatives**. Presence was benchmarked; a fully fitted two-part model was not.

## 7. Biological findings and confidence audits

Four descriptive WT neighborhoods (20-bp summaries around smoothed peaks) have small-call coverage 6.51%, 5.05%, 4.95% and 2.26%: **P1 chr19:11089316–11089335; P2 11089338–11089357; P3 11089378–11089397; P4 11089413–11089432**. These are data-selected neighborhoods, not inferred protein contact boundaries.

**Position 125 / chr19:11089389:** all substitutions reduce full-WT P3 coverage from 4.95% to **1.37%, 1.23%, 1.39%** (C/G/T; 380/480/438 barcodes). Effects persist across barcode halves, equal molecule weighting, depth restrictions, leave-one-barcode influence checks, distal methylation/protection standardization, and shared A/T positions excluding mutation ±2 bp. Those initial shared-site analyses reused caller annotations.

The latest audit independently examines **retained m6A calls** at nine unchanged P3 A/T positions outside ±2 bp. Barcode-weighted changes are **+4.12/+3.65/+7.00 pp** for C/G/T; both genomic adenine strands increase. Whole-element raw changes are smaller, +0.93/+0.34/+2.33 pp. This supports a localized change, including for A>T, which preserves A/T count but switches strand. It does not eliminate sequence-context or caller-boundary artifacts.

Barcode-independent direct-genotype comparisons retain the pattern. All 50 one-molecule-per-barcode selections retain positive raw and negative footprint changes. ≥40-pass filtering attenuates the footprint loss and leaves only 158–234 molecules per variant; it is a limitation, not a filter to optimize away. Published activity is not uniformly reduced (+0.08/+0.47/−0.07 in source 1), so this is **not an identified activating TF footprint**.

**Adjacent region / chr19:11089407–11089414, `CTCCTCCC`:** in the full-WT atlas, 22/24 SNVs reduce promoter MSP (median −2.00 pp), 22/24 increase 85–250-bp protection (+2.28 pp), and 23/24 reduce activity in each published source (median log2 −2.50/−3.08). Local TF changes are less coherent, with only two A/T opportunities in eight bases. The latest raw audit finds lower whole-element retained methylation in 21/24 and lower MSP in 20/24; barcode-independent grouping gives 23/24 and 22/24. The group trend survives; individual near-zero signs are less stable.

**Do not conflate WT references:** the full atlas uses 268,283 WT barcodes. The raw audit uses the fixed **10,000-barcode /37,075-molecule** WT sample, or a separate deterministic 1/32 read sample for genotype-only comparisons. This explains small changes in reported differences and adjacent-region direction counts. All qualifying candidate molecules are retained. In the local 143–150 interval, masking ±2 bp leaves **zero raw opportunities for variants at 145/146**; those endpoints are missing.

The latest 18 settings include weighting, barcode independence, circularity/length, read quality/passes/mapQ, alignment orientation, barcode depth, baseQ, retained ML≥250 and whole/distal methylation restrictions. MapQ≥60 and whole m6A<.90 retain the original cohort and add no evidence. Noncircular subsets contain only 18–93 molecules per candidate. Distal conditioning can remove genuine biological effects. ML/baseQ sensitivities do not re-call footprints. WT subsampling ranges are descriptive sampling scales, **not confidence intervals or p-values**.

**Motifs:** 879 JASPAR2024 and 1,443 HOCOMOCOv12 matrices were scanned locally on both strands, with rounded log-odds dynamic-program tails and an unadjusted 1e-4 sequence threshold. There are 114 matches across 923,098 windows (110 pass GC-background sensitivity), not 114 binding discoveries. No variant-created-site scan was done. Fixed-site ALT scores cover 5,022 motif/variant pairs.

The known SRE1 `ATCACCCCAC` maps to **chr19:11089391–11089400**, two bases downstream of position 125. SREBF motif-score changes correlate with published activity (Spearman .76/.72 for SREBF1), but hardly with local small-call coverage (.02). An overlapping SP1 model correlates with MSP (.50), activity (.69/.75), and less with local TF (.21). Overlapping SP/KLF/WT1 models prevent unique attribution. These post hoc associations motivate experiments; they do not identify binding proteins, causal mediation, activators or repressors.

## 8. Reproduce, validate or extend

The [portable evidence index](evidence/20261020/README.md) includes selected configurations, audit/release records, calibration summaries, coordinate evidence, candidate tables and checksums. It is an inspection package, not a substitute for raw data or all large intermediates. `results/` is Git-ignored. Full local outputs are preserved under the exact dated paths in the reports.

Key entrypoints, in dependency order:

| Stage | Workflow / code | Instructions |
|---|---|---|
| Original pilot catalog / published function joins | `workflow/pilot.smk`, `fsmpra.catalog` | [Pilot usage](pilot_usage.md) |
| New raw BAM through features | `workflow/process.smk`, `config/ldlr_processing.json` | [Processing](processing_workflow.md); only bounded smoke test completed end-to-end |
| Existing-alignment complete ledger | `workflow/ledger.smk`, `fsmpra.ledger` | [Ledger usage](rep3_ledger.md); saved run config in evidence bundle |
| Reuse historical calls | `workflow/footprint_join.smk` | [Join usage](footprint_join.md) |
| Barcode arrays, diagnostic inference/calibration | `workflow/statistics.smk` | [Statistics usage](statistics_usage.md) |
| Molecule heterogeneity and barcode summaries | `workflow/heterogeneity.smk` | [Heterogeneity usage](heterogeneity_usage.md) |
| WT sampled raw calls, quality, first benchmark | `workflow/footprint_quality.smk` | [Quality usage](footprint_quality_usage.md) |
| Large WT presence/coverage and calibration | `workflow/heldout_calibration.smk` | [Held-out usage](heldout_calibration_usage.md) |
| Landscape and genomic atlas | `workflow/landscape.smk`, `fsmpra.genomic_coordinates` | [Landscape](landscape_usage.md), [coordinate reproduction](genomic_coordinate_audit_20261020.md) |
| Motifs and barcode audits | `fsmpra.biology`, `tools.render_biology_atlas` | [Candidate report](rep3_candidate_biology_20261020.md), [viewer usage](biology_atlas_usage.md) |
| Candidate raw extraction/sensitivities/viewer | `fsmpra.raw_candidate_extract`, `fsmpra.raw_candidate_summary`, `tools.render_raw_candidates` | [Raw-audit commands](rep3_raw_candidates_20261020.md) |

The latest complete viewers are `results/rep3_landscape_20261018/atlas_biology_20261020/index.html` and `results/rep3_raw_candidates_20261020/viewer/index.html`. They embed data and require no server; linked exports require the accompanying directory. The biology source directory is `biology_20261020_complete`, not its earlier development sibling. Viewer templates and generators are in Git; generated HTML is local.

Use `pixi.toml`/`pixi.lock` as environment references, but historical stages ran existing local Python/Snakemake/tool environments, not one newly verified portable environment. Tool binaries/models and project data paths are installation-specific. Configs retain exact references; the main BAMs are the `PS01898.hg38_LDLR.ft.aligned.bam` alignment and `plasmid_fs/footprints/PS01898_LDLR_HepG2.fiberhmm.bam` under the rep3 data directory. The design file is under the lab's `APC_Promoter_SGE/probe_design/LDLR_oligos.txt`. Local Python was `/mmfs1/home/bmallo/miniconda3/bin/python`; Snakemake was in its `envs/snakemake` environment.

Run substantive work only inside a verified Slurm allocation. Default these to one and request resources through Slurm, never from visible node totals:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export MPLCONFIGDIR=/tmp/fsmpra-matplotlib
python -m unittest discover -s tests -p 'test_*.py'
```

Ledger extraction used four allocated CPUs; most subsequent stages used one. Raw genome/plasmid alignment needs substantially more RAM than downstream summaries. Follow per-stage resource declarations; do not blindly invoke the legacy Pixi eight-core task. Use node-local scratch for temporary SQLite/BAM work. No server is needed for the viewers.

**Reproduction pitfalls:** outputs are atomic and refuse overwrite; use new directories and copied configs. Later stages reference dated upstream paths, and some candidate tools are rep3-specific. Changing `--output` alone does not retarget every dependency. Frozen source guards check hashes and/or path/size/mtime; a clean Git checkout changes timestamps, so historical audits may fail despite identical file contents. Recreate audits from verified inputs or reproduce from preserved snapshots—do not silently disable guards or edit historical audits to force a pass. Large BAM provenance is generally metadata, not full-content hashing. Some stages used direct CLIs with subsequent Snakemake dry runs; audit/log records, not absent Snakemake execution metadata, document those runs.

Independent checks should prioritize read-ID/sequence/coordinate reconciliation, per-barcode depths, missing-vs-empty behavior, and recomputation of original TF/MSP means before new inference. Existing validators are `tools/validate_*`, with focused unit tests under `tests/`. Full heterogeneity and held-out endpoint checks matched prior barcode arrays to 2.22e-16; landscape reconciliation covered all prior profiles. The latest raw audit validates 15,552 comparisons, 106,344 per-base rows, 504 retention rows, joint-state mass and original variant profiles. Nine focused raw tests and all 1,944 variant×setting×region DOM/canvas checks passed. Scientific PNGs were inspected; **full browser visual QA was not performed**. Software agreement does not establish biological validity.

## 9. Useful next work that has not already been settled

1. **Focused experimental replication:** use the proposed ten-construct panel (WT, all 125 substitutions, GC-preserving adjacent variants, A/T strand swaps and 147T>G functional exception), balanced barcodes and independent preparations. Add separately prepared sequence-matched protein-free reactions and a documented labeling-time/enzyme series. These controls have not been generated.
2. **Measurement calibration:** recover original kinetics/probabilities for new data; test context/strand-specific labeling and calling on known inputs. Compare alternative caller settings or caller models on matched molecules, keeping total and conditional effects separate. Full-data recaller benchmarking and sequence-matched controls have not been completed.
3. **Targeted statistical development, if needed:** predefine a small candidate family on independent data; evaluate fresh permutation banks, rare-event estimators and calibrated hierarchical/two-part models or interval inversion. Benchmark unequal dispersion, realistic dependence and gains/losses. Do not re-present Welch, ordinary barcode permutation, presence outcomes, distal stratification or a one-read-per-barcode audit as untried ideas.
4. **Functional integration:** import new MPRAsnakeflow barcode-level RNA/DNA counts, metadata/QC and activity estimates with library/sequence identity validation. Current published joins do not provide matched-sample mediation or protein regulatory-sign inference.
5. **Motif follow-up:** independent HepG2 expression/binding evidence, candidate-specific perturbations or orthogonal assays; variant-created-site scanning or response-pattern motif inference remain possibilities. Existing work scanned WT reference motifs and fixed-site substitution scores, not de novo binding motifs inferred from response clusters.
6. **Engineering generalization:** fresh-sample end-to-end execution, portable path/config management, integrated Snakemake rules for the latest biology/raw-audit tools, automated source portability, and full browser QA. The present scripts are auditable research stages, not a finished generic production/dashboard service.

Do not pool different libraries/reaction conditions, label high-methylation molecules as known naked controls, equate long annotations to single nucleosomes, convert unavailable measurements to zero, tune filters against functional agreement/significance counts, or interpret the withheld candidate q-values as discoveries.
