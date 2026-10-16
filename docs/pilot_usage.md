# Executable redesign pilot

This is the first implementation of the redesign. It validates the library and
functional-data joins, imports an auditable molecule sample, and measures
small-footprint occupancy and an accessibility proxy. It does **not** yet make
statistical discovery calls or replace the legacy workflow.

## Run

Run from the repository root, inside a Slurm allocation for the molecule and
BAM scans. Python 3.9+ and `pysam` are sufficient; the catalog and molecule
import use only the standard library. No environment installation is required
when those are already available. Every command uses one processing thread;
the BAM reader explicitly uses one thread.

```bash
python -m fsmpra catalog \
  --manifest config/ldlr_pilot.json --output results/pilot/catalog

python -m fsmpra molecules \
  --catalog results/pilot/catalog/catalog.json --sample HepG2_rep3 \
  --fraction 0.01 --seed 42 --max-molecules 200000 \
  --output results/pilot/HepG2_rep3/molecules

python -m fsmpra measure \
  --catalog results/pilot/catalog/catalog.json \
  --molecules results/pilot/HepG2_rep3/molecules \
  --output results/pilot/HepG2_rep3/measurements

python -m fsmpra dashboard \
  --pilot results/pilot/HepG2_rep3/measurements/pilot.json \
  --output results/pilot/HepG2_rep3/dashboard

# Optional static overview (requires numpy and matplotlib):
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m fsmpra plot \
  --pilot results/pilot/HepG2_rep3/measurements/pilot.json \
  --output results/pilot/HepG2_rep3/plots
```

Outputs cannot overwrite an existing directory. Use a new run directory when
inputs or choices change. Source data are read-only. Each stage checks its
input provenance; results include source paths and versions. Small inputs are
SHA256-hashed; large-source identity uses path, size, and modification time.

An independent Snakemake entry point is also provided:

```bash
snakemake -s workflow/pilot.smk --cores 1 --config \
  pilot_run_dir=results/pilot_run pilot_sample=HepG2_rep3 pilot_fraction=0.01
```

Use the site's current Slurm configuration to select account/partition/memory
for job submission. No historical partition or account is hard-coded here.
The existing `workflow/Snakefile` continues to run the old analysis.
If Snakemake lives in a different environment from pysam, pass
`pilot_python=/absolute/path/to/python-with-pysam` through `--config`.

## What is imported

`config/ldlr_pilot.json` explicitly names the newer HepG2 library, its pinned
reference sequence, design oligos, the three independent biological replicates,
and the two external published MPRA tables. The HEK293T library is not assigned
the newer library's coordinates. Cross-library comparisons remain pending.

The catalog reconstructs the sequence of every designed WT/SNV/one-base indel,
collapses sequence-equivalent oligos, and assigns a library-qualified sequence
hash as the stable ID. Original oligo names, canonical edit labels, and
equivalent indel positions remain available. SNV positions are stored explicitly
as zero-based rotated-reference positions and one-based oligo/genomic positions.
An insertion position denotes its right-flank boundary, not a substituted base.

MPRA records retain source identity, original effects and p-values, and join
status. Equivalent deletion measurements remain separate records; they are not
averaged. The source tables' aggregate DNA/RNA columns are preserved, without
re-estimating activity or treating them as biological replicate count matrices.

Historical coverage is imported with `historical_exploratory` status and its
original cohort names. `final` denotes the old analysis's selected cohort; it
is not the new pipeline's recommended filter.

## Molecule sample

Read selection uses a stable 64-bit BLAKE2b hash of seed and read ID across the
**entire genotype ledger**. Selection is independent of file order, methylation,
or genotype, and is nested when the fraction increases with the same seed.
It is a Bernoulli-style fraction, not a requested exact number of reads. These
are sampled plasmid molecules, not a random sample of all delivered reads.

The importer scans the complete compressed ledger and metrics files. It retains
all selected genotype classes, missing metric rows, barcodes, source confidence
tiers, global methylation, and original pipeline filter annotations. Duplicate
selected keys and stale sources cause errors. No methylation filter is applied.
An explicit selected-read cap bounds memory; exceeding it fails rather than
silently truncating or changing the sampling scheme. The pilot is intended for
bounded samples; full-library storage needs a subsequent disk-backed adapter.

## Initial measurements

The user-confirmed primary analysis interval is the **211-bp mutable element**.
Default measurement cohort: upstream
`plasmid_full`, WT or designed single SNV, at least 10 polymerase passes, and
exact agreement of the BAM's mutable sequence with the assigned design. All
methylation levels remain included. Those quality choices are configurable and
exploratory; the import preserves the data needed to change them.

Measurements use the existing FiberHMM `MA` tracks in stored-SEQ coordinates:

- `tf_small`: occupancy of called TF-track segments 10–60 bp (configurable).
- `nucleosome`: occupancy of nuc-track segments at least 85 bp (configurable).
- `msp`: coverage by caller-defined accessible segments, used as an initial
  accessibility **proxy**.

Each position has separate called-read and eligible-read counts. A missing
`MA` annotation is missing data. FiberHMM's sparse writer omits empty track
sections: for a sample explicitly declaring the evaluated tracks in its
`footprint_contract`, an omitted evaluated track contributes zero. Without
that declaration, an absent track remains missing, and the audit counts it.
An explicitly present empty track always contributes measured zero. Query
intervals project only onto aligned bases, without bridging reference deletions.
Length distributions preserve the original query-space footprint lengths.

Per-molecule MSP fractions and histograms distinguish changes in mean coverage
from changes in the molecule distribution. Raw WT/variant profiles are retained.
An additional common-site summary excludes the variant position ±2 bp from
both the variant and its WT comparator; this is a sensitivity readout, **not a
validated correction for sequence-dependent detectability**. Small-footprint
and accessibility measurements are reported independently of MPRA activity.

The pilot emits no p-values, confidence intervals, inferred TF identities, or
claims of significant chromatin effects. A/T-dependent observability, caller
calibration, protocol differences, biological replication, and the appropriate
statistical family still require evaluation. Runs 1–2 used 10 minutes and
1.5 µL Hia5 per 2 million cells; run 3 used 5 minutes and 3 µL per 2 million
cells. They must not be pooled as three exchangeable biological replicates.

## Compare reaction conditions

After producing one measurement bundle per sample, run:

```bash
python -m fsmpra compare --manifest config/ldlr_pilot.json \
  --pilots results/conditions_20261014/HepG2_rep1/measurements/pilot.json \
           results/conditions_20261014/HepG2_rep2/measurements/pilot.json \
           results/pilot_20261014/HepG2_rep3/measurements/pilot.json \
  --min-reads 10 --output results/my_comparison

python -m fsmpra dashboard --pilot results/my_comparison/collection.json \
  --output results/my_comparison_dashboard
```

Alternatively edit `config/ldlr_comparison.json` and run
`snakemake -s workflow/compare.smk --cores 1`; set `pilot_python` through
`--config` if Snakemake and pysam use different environments.

Comparison rejects differing references, intervals, allele sets, measurement
parameters, missingness conventions, repeated sample IDs, and repeated BAM
sources. Updated experimental annotations are applied to a new collection,
with the annotation manifest recorded in provenance; original results remain
unchanged. The comparison verifies source identities before reading them.

Outputs are `collection.json` (dashboard data), `comparison.json` (QC,
processing evidence, and descriptive pairwise summaries), `sample_qc.tsv`,
and `allele_comparison.tsv`. Conditions remain separate. Pairwise correlations
use only SNVs meeting the minimum coverage in both runs; absent/constant or
fewer-than-three valid pairs produce null. No significance test is performed.

The comparison dashboard includes reaction labels, coverage and full-plasmid
fractions, WT profile overlays, and a selected-allele table across runs. The
full-plasmid fraction's denominator is the metrics-matched **sampled genotype
ledger**, not the total sequencing output. Low-depth measurements and differences
in reaction conditions/software prevent a causal interpretation of between-run
differences.

## Output contract: fsmpra.pilot.v1

| Artifact | Purpose |
|---|---|
| `catalog/catalog.json` | References, designed sequences, sample provenance, and published MPRA records |
| `catalog/variants.tsv` | Stable variant IDs and explicit coordinate/allele fields |
| `catalog/mpra.tsv` | External measurements and join status, one row per source record |
| `catalog/historical_coverage.tsv` | Existing coverage with original filter/cohort labels |
| `molecules/molecules.tsv.gz` | Selected ledger rows joined to source metrics; missing values are empty |
| `molecules/audit.json` | Whole-ledger class counts, sample fraction/seed, joins, and provenance |
| `measurements/molecule_features.tsv.gz` | Per-read footprint/MSP fractions and global methylation |
| `measurements/variant_summary.tsv` | Separate WT/SNV measurements and descriptive differences |
| `measurements/pilot.json` | Dashboard-ready profiles, counts, distributions, sampled molecules, and MPRA context |
| `measurements/audit.json` | Cohort exclusions, parameters, source identity, and limitations |
| `dashboard/index.html` | Portable interactive viewer with embedded data, styles, and scripts |

JSON uses `null` for unavailable values and rejects NaN/Infinity. Every
designed SNV appears in the result, including those with no sampled evidence.
Molecule examples use deterministic bottom-k read hashes, not the first reads
in the BAM. Example intervals are zero-based, half-open, relative to the
declared analysis interval.

## Local dashboard

Open the exported `dashboard/index.html` in a browser; it needs no server,
installation, network connection, or data upload. The overview links alternate
alleles to WT/variant footprint profiles, molecule-level accessibility and
methylation distributions, deterministic molecule examples, and separate
published MPRA records. The coverage threshold changes the overview display,
not the measurement cohort. Differences remain descriptive.

Use **Load pilot JSON** to inspect additional measurement bundles locally.
Replicates remain separate datasets. Download a selected allele's evidence,
the summary CSV, or the displayed profiles as SVG. The measurement definitions
and source provenance travel with the HTML; it contains read identities and
local source paths, so treat it as a research data artifact when sharing.

The source template in `dashboard/` is bundled by the CLI; use the exported
HTML to view a measured dataset. A production comparison interface, calibrated
statistics, and cross-replicate analysis remain later stages.

## Tests

```bash
python -m unittest discover -s tests -v
node tests/test_dashboard.js
```

The synthetic end-to-end tests exercise a reverse-strand WT, a SNV, methylation
extremes, missing metrics, empty versus absent footprint tracks, coordinate
projection, and allele-aware MPRA joins. Additional tests check equivalent
indels, input disagreements, duplicate joins, and sampling limits.
