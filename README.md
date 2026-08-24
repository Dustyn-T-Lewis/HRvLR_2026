# A_HRvLR_2026_V2

Maps the skeletal-muscle proteome onto training adaptation continuously. No
responder groups, no cut points, no baseline contrast.

The original study labelled subjects High and Low Responder by median-splitting
a composite hypertrophy score. V1 tested that label at four levels and found
nothing. This project drops the label and asks the question the label was
standing in for: **does a protein's change over a window track how much a
subject adapted?**

## Layout

```
00_input/          raw matrix, sample metadata, phenotype table
01_Filtering/      protein filtering
02_Normalization/  cycloess, imputation, per-sample pathway scores
03_Features/       modules, the association sweep, its calibration
04_Figures/        F01_phenotype, F02_association
functions/         shared code, flat
```

Every stage carries `a_script/` `b_reports/` `c_data/`; pipeline steps are
numbered `NN_*.R` and figure units use `setup.R` + `01_run_<name>.R` +
`composite.R` + `panels/`. `c_data/` is tracked, `b_reports/` renders are not,
paths go through `here::here()`, and any stochastic step is seeded.

## The design

Each subject contributes one column: how much a feature moved over a window.
That change is regressed on the subject's adaptation. Nobody is cut into a
group, so no cut point has to be defended and no composite can separate its own
ingredients.

- **Windows.** Training (T2 − T1) and Acute (T3 − T2). Baseline is absent by
  design: comparing levels between people answers a different question from
  whether a change tracks a change, and V1 already tested the baseline form
  across 54 cells without promoting anything.
- **Levels.** Proteins (1900), WGCNA module eigengenes (12), singscore
  pathways (57).
- **Estimator.** `limma` with a continuous predictor. One row per subject means
  no repeated measures inside the fit, so no blocking and no
  `duplicateCorrelation` — the within-subject structure is spent forming the
  difference. The moderated variance is why limma beats a per-feature `lm` at
  n = 14.

3 levels × 2 windows × 10 phenotypes = **60 cells**. BH within each cell, never
across them.

## Phenotypes

Ten, up from the five V1 used. The additions were already in the input and no
stage had read them.

- `comp_hypertrophy`, the composite the original label was cut from.
- `d_fcsa_I`, `d_fcsa_II`, `d_fcsa_mixed`, fibre cross-sectional area.
- `d_nfibre_mixed`, `d_nfibre_I`, the MyoVision columns. These are **fibre
  counts, not areas**, despite carrying "fCSA" in their meta names; the source
  workbook calls them "Number of fCSA - Mixed (MyoVision)". They run 142-1119
  where the areas run 3800-10700 and move against area, because larger fibres
  pack fewer into the imaged field.
- `d_mcsa`, whole-muscle cross-sectional area.
- `d_1rm_legpress`, `d_1rm_ext`, strength.
- `volume_load`, total kilograms lifted, one value per subject, spanning
  3.5-fold. The only variable describing what a subject did rather than what
  happened to them.

They are not ten independent questions. Three fibre-area measures correlate
above 0.9 with each other and 0.85 with the composite; the two fibre counts
correlate 0.94 with each other and −0.5 to −0.7 with area. That leaves roughly
five independent axes: fibre size, whole-muscle CSA, each 1RM, and volume load.
F01 panel B is that structure.

Whole-muscle CSA and both strength measures rose over training (d = 1.13 to
1.48). No fibre-area or fibre-count measure moved at all, every interval
covering zero. Phenotype exists at T1 and T2 only, so with no comparator arm
and no repeat baseline, true individual response cannot be separated from
measurement error within this study (Atkinson & Batterham 2015, *Exp Physiol*
100:577).

## What it found

**One cell of 60 cleared BH, against 2.6 expected by chance.**

| Level | Cells | Cleared BH | Expected under the null |
|---|---|---|---|
| Proteins | 20 | 0 | 0.77 |
| Modules | 20 | 1 | 0.92 |
| Pathways | 20 | 0 | 0.92 |
| **All** | **60** | **1** | **2.61** |

The expectation is not one rate applied sixty times. BH over 12 module
eigengenes is a far weaker filter than the same alpha over 1900 proteins, so
each cell's null rate is estimated separately by shuffling the phenotype across
subjects, 999 times per cell.

**The one hit does not hold.** Module greenyellow against Δ whole-muscle CSA
over the acute window: BH = 0.037, n = 15. Three independent checks disagree
with it.

- Its own permutation null puts it at p = 0.051.
- Spearman rho = −0.35, p = 0.20. A linear fit a rank test cannot see is being
  carried by the extremes of the scale rather than by the ordering.
- Refit dropping each subject in turn, it holds in **6 of 15 folds**. Dropping
  HR_S29 or LR_S14 moves it to BH ≈ 0.20.

LR_S14 is the same subject V1's `06_mcsa_axis` found its Δ mCSA result resting
on. Two independent analyses, one influential point.

**V1's one hit does not reproduce, and could not have.** PSME1's best cell here
is p = 0.035 unadjusted, BH = 0.885. V1 tested its T2 *level* against Δ mCSA;
this tests its *change*. Different quantities, so this is not a failed
replication.

## Reading the null

Nothing in this proteome tracks how much these subjects adapted, at any of
three feature levels, over either window, against any of ten phenotypes. The
sweep returned fewer hits than chance predicts, and the single hit fails
permutation, rank correlation, and leave-one-subject-out.

That is a calibrated negative rather than an absence of evidence. Each stage
carries its own null, so "we found nothing" and "there was nothing to find" are
distinguishable here.

## Running it

```r
source("setup.R")
for (f in list.files("03_Features/a_script", "[.]R$", full.names = TRUE)) {
  source(f)
}
for (f in list.files("04_Figures", "^01_run_.*[.]R$",
                     recursive = TRUE, full.names = TRUE)) {
  source(f)
}
```

`03_confirm_hits.R` permutes all 60 cells and takes about 45 minutes.
Everything else runs in under a minute.

## Archived

`archive/` holds the classification work this project started as: a sweep of
ten candidate group labels across three contrasts, an unsupervised subtype
search, and their figures. It is untracked and superseded. Its conclusions ran
the same way — the label separates the one measure that never changed, the
proteome carries no subtype structure that beats a no-cluster null, and 2 of 24
label cells cleared BH against 2.6 expected.

`fgsea` was dropped with it. On this proteome a random relabelling of the same
subjects produced a median of 102 significant sets against 98 observed
(p = 0.52), because its preranked null permutes gene labels and so treats
co-regulated proteins as exchangeable. Pathway work here goes through
singscore, which scores each sample independently and is fitted through the
same estimator the proteins use.
