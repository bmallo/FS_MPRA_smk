# Replicate-3 statistics stage

The fixed protocol is in [rep3_statistics_design.md](rep3_statistics_design.md).
Run from the repository root inside a Slurm allocation, using one CPU per rule:

```bash
snakemake --snakefile workflow/statistics.smk --cores 1 \
  --resources mem_mb=8192
```

The default output is `results/rep3_statistics_20261017/`. Each stage directory
is published atomically and refuses to replace existing results. `aggregates/`
contains `barcode_counts.npy` with shape `(N,2,211)` and uint32 base counts in
`tf_small,msp` order; `barcode_molecules.npy` and `barcode_metadata.tsv.gz`
identify each row. `cohort_counts.npy` and `cohort_molecules.npy` are catalog-order
uint64 sums for all eligible molecules, including unsupported barcodes.
`variant_metadata.tsv` fixes their catalog-order IDs. The join is read once.

`inference/effects.tsv.gz` has 13,926 fixed SNV × endpoint rows. Effects are
variant minus WT fractions, with percentage-point columns. Testable rows have
candidate Welch p-values and joint BY q-values. Untestable rows retain means
and status but have blank uncertainty and candidate p/q fields. The family
still counts each such member internally as p=1. `profiles.tsv.gz` retains
211-position TF/MSP descriptive curves. `sensitivities.tsv.gz` records
supported-cohort molecule weighting, all-eligible molecule weighting, and
mutation ±2 masking. Masked results omit the same reference bases in both
arms. No profile, masked value, or molecule-weighted value inherits the fixed
window family error control.

`release_status.json` leaves statistical support false until a calibration
report with `status: "pass"` and `validated_endpoint_ids` is supplied to
`python -m fsmpra.statistics infer --calibration REPORT.json` for a new output
directory. A row's `caller_effect_supported` requires a validated endpoint,
testable contrast, and candidate q ≤ 0.05. Even after statistical calibration,
sequence-aware caller detectability and biological replication remain
unresolved. Candidate p/q values describe approximate within-sample caller
associations, not validated binding changes.

The workflow runs calibration between aggregation and inference by default.
Set `--config statistics_run_calibration=false` to produce candidate inference
without a calibration report; release status then records `absent` and withholds
supported effects.

`calibration/report.json` records all 88 WT A/A strata (four sample sizes ×
22 endpoints) and 27 simulation scenarios, with 200 iterations each.
`calibration/iterations.tsv.gz` preserves individual tests. Any assessed
failure blocks release. Unresolved null simulation scenarios also block a
global pass until inference can be restricted to a validated outcome regime.
The exact binomial intervals measure Monte Carlo uncertainty under these
checks; they are not confidence in biological replication or sequence bias.

The first endpoint measures mean fraction of a window covered by 10–60-bp
caller intervals. It combines changes in call prevalence, position and length;
it is not yet a binary occupancy estimate for a particular TF-binding site.
The promoter endpoint is mean MSP-covered fraction over the 211-bp element.
Both remain summaries of the supplied caller output.

To preserve existing run directories, choose a fresh output for a new analysis:

```bash
snakemake --snakefile workflow/statistics.smk --cores 1 \
  --resources mem_mb=8192 \
  --config statistics_output=results/my_new_statistics_run
```

The complete manifest, code/config snapshots and logs for the first execution
are in `results/rep3_statistics_20261017_run/`. The joined table is unchanged.
