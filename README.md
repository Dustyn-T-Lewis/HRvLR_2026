# A_HRvLR_2026_V2

A second pass at the HR/LR skeletal-muscle proteomics dataset, built on V1's
preprocessing and asking two questions V1 never ran: what the responder label
is made of, and whether the proteome groups these subjects on its own.

Stages 00 through 02 are copied from V1 unchanged. `03_Features` runs in
dependency order; each stage reads the previous stage's `c_data/` and writes
its own.

```
00_input/  01_Filtering/  02_Normalization/     copied from V1, not re-derived
functions/                                      shared code, flat
03_Features/
  01_Responsiveness/   label audit, candidate label construction
  02_WGCNA/            modules and eigengenes
  03_Subtypes/         blind clustering with a no-cluster null
  04_Proteins/         limma sweep across candidate labels, then permutation
  05_Pathways/         fgsea and singscore, then a calibration check
04_Figures/
  F01_phenotype/  F02_subtypes/  F03_proteome/
tests/
```

Conventions follow V1: `a_script/` `b_reports/` `c_data/` at every unit,
numbered `NN_*.R` for pipeline stages, `setup.R` + `01_run_<name>.R` +
`composite.R` + `panels/` for figure units, one source-data workbook per
figure, `c_data/` tracked and `b_reports/` renders ignored, `here::here()` for
every path, seed 42 before any stochastic step.

## Declared gates

Written before fitting. A stage that fails its gate reports the failure and
stops rather than reinterpreting.

1. Stage 03 interprets clusters only if BIC selects more than one component or
   the observed fit beats a no-cluster null. **Shut.**
2. Stage 05 runs only if stage 04 produces a BH survivor. **Open** (three
   survivors), so stage 05 ran.
3. Permutation runs only to confirm a hit that already exists, never as a
   standing sweep. **Armed twice**, at stage 04 and again at stage 05.

## What it found

**The label.** HR/LR is the exact median cut of `comp_hypertrophy`: the top
eight subjects are HR without exception. That composite arrives without a
stated formula but is reconstructable from the five measured outcomes at
r² = 0.977, and it is dominated by the two fibre cross-sectional area measures
(r = 0.86 and 0.83) over whole-muscle CSA (0.50) and the two strength measures
(0.27, 0.13).

Neither fibre measure changed over training: d = −0.11 [−0.64, 0.43] for type
I and 0.06 [−0.47, 0.60] for type II. Whole-muscle CSA and both strength
measures did, at d = 1.13 to 1.48. The label separates the fibre measures
(p ≈ 0.001) and, weakly, whole-muscle CSA (p = 0.044); it does not separate
either strength measure (p = 0.79 and 0.63). Splitting on a composite
necessarily separates that composite's ingredients, so the fibre separation
says nothing about whether the label tracks anything beyond its construction.

The composite is weakly bimodal — bootstrap likelihood-ratio test over 999
replicates gives p = 0.035 — and that two-component solution recovers the given
label exactly.

**Subtypes.** No. Model-based clustering of the baseline proteome, blind to
every label, was run on module eigengenes and on the 500 highest-variance
proteins, at two, three and four retained principal components. No cell beat
its null (p = 0.10 to 0.61). The null matters: data drawn from a single
Gaussian with the observed covariance and no groups in it gets assigned more
than one component in 47 to 77 percent of draws at this sample size. A forced
two-group split agrees with no label (adjusted Rand index −0.08 to 0.16).

**Proteins.** Six labels by two contrasts is twelve cells. Three
protein-contrast survivors at BH < 0.05, all at Baseline: TMED5 and TBCB under
whole-muscle CSA, PLIN3 under 1RM leg extension. Neither hit cell survives a
999-permutation subject-label null (p = 0.056 and 0.130) — a random split of
these sixteen subjects clears BH 13 to 16 percent of the time. The given HR/LR
label is the most null of the six, at a smallest adjusted p of 0.996 and 0.938.

**Pathways, and a warning.** fgsea returned 1,398 set-contrast hits down to
padj = 2e-19. It returned more under 1RM leg extension, a split with almost no
agreement with the responder label, than under the responder label itself. A
calibration run says why: across 100 random relabellings of the same subjects,
the median number of significant sets is 102, against 98 observed
(empirical p = 0.52). fgsea's preranked null permutes gene labels and so treats
proteins as exchangeable; pathway members are co-regulated and share technical
structure in a normalised MS matrix, which makes that null anticonservative
here. **No pathway result in this repository should be read as evidence.**
singscore, scored per sample and fitted through the same estimator the proteins
used, returned zero hits across all twelve cells. V1 recorded the same failure
mode for STRING's PPI enrichment p on this proteome.

## Equivalence with V1

The sweep estimator generalises V1's design to an arbitrary label. Under the
given HR/LR label its two contrasts are V1's `Baseline_HRvLR` and
`Training_Interaction`, and they agree with V1's committed numbers to
4.9e-15 across all 1,900 proteins. `verify_v1_equivalence()` runs this check
and its result ships in `03_Features/04_Proteins/c_data/01_sweep.xlsx`. Module
eigengenes reproduce V1's to zero difference, as they must — module
construction never sees a label.

## Data integrity note

Manual and MyoVision fibre cross-sectional area correlate negatively
(r = −0.63 for mixed, −0.44 for type I) with a mean offset of −5,660 across 32
paired measurements. Two methods measuring one quantity cannot correlate
negatively, so either they measure different things despite the naming or
something is mislabelled upstream. Nothing here uses the MyoVision columns —
`build_phenotype.R` never did — but the discrepancy is recorded in
`03_Features/01_Responsiveness/c_data/01_label_audit.xlsx` because it touches
the measure the responder label is built on.

Phenotype exists at T1 and T2 only. With no comparator arm and no repeat
baseline, true individual response cannot be separated from measurement error
within this study (Atkinson & Batterham 2015, *Exp Physiol* 100:577).

## Deliberately not built

Supervised screens. V1 ran elastic net with nested leave-one-subject-out
across 12 classification cells and 54 continuous association cells and promoted
nothing.

A new latent phenotype model. V1's GALAMM pilot converged at n = 16 and found
the latent factor is `comp_hypertrophy` itself (rho = 0.94), with zero protein
survivors at a smallest q of 0.117.

`d_mcsa` as an external validator. V1 established it is an ingredient of the
arm definition, roughly 36 percent of the composite's formula.

## Running it

```r
source("setup.R")
for (f in list.files("03_Features", "^[0-9]{2}_.*[.]R$",
                     recursive = TRUE, full.names = TRUE)) source(f)
for (f in list.files("04_Figures", "^01_run_.*[.]R$",
                     recursive = TRUE, full.names = TRUE)) source(f)
```

The two permutation stages (`04_Proteins/02_confirm_hits.R`,
`05_Pathways/02_fgsea_calibration.R`) take several minutes each. Everything
else runs in under a minute.
