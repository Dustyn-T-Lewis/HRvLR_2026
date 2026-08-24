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

- **Windows.** Six. Three levels (the value at T1, T2 or T3) and three changes
  (training T2−T1, acute T3−T2, total T3−T1). Levels ask a between-person
  question, changes a within-person one; the two families are reported apart
  and never pooled.
- **Levels.** Proteins (1900), WGCNA module eigengenes (12), singscore
  pathways (57).
- **Estimator.** `limma` with a continuous predictor. One row per subject means
  no repeated measures inside the fit, so no blocking and no
  `duplicateCorrelation` — the within-subject structure is spent forming the
  difference. The moderated variance is why limma beats a per-feature `lm` at
  n = 14.

3 feature levels × 6 windows × 10 phenotypes = **180 cells**. BH within each
cell, then a second correction across all 180.

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

**Nothing survives the sweep-level correction.** Three of 180 cells clear BH
inside themselves, all against change in whole-muscle CSA. After correcting the
per-cell permutation p across the sweep, none has q below 1.

Every count runs below its own null:

| | observed | expected under the null |
|---|---|---|
| Cells with permutation p < 0.05 | 3 | 9 |
| Cells with ≥1 BH hit | 3 | 8.9 |
| Total BH hits | 6 | 24.4 |
| Surviving BH across the 180 cells | **0** | — |

The expectation is built per cell, by shuffling the phenotype across subjects
999 times. It has to be: BH over 12 module eigengenes is a far weaker filter
than the same alpha over 1900 proteins, so a single rate applied 180 times
would be wrong in both directions.

### The one that came closest

Module **greenyellow** against change in whole-muscle CSA, at the T2 level.
It passes every check applied to itself and still fails the sweep.

- BH = 0.0009 within its cell; **0.0058** after adjusting for the biopsy
  composition panels that confound the phenotype.
- Spearman rho = 0.59, p = 0.021, so a rank test sees it.
- Refit dropping each subject in turn, it holds in **14 of 15** folds.
- Its own permutation p is 0.020 — but **q = 1** across the sweep, and three
  cells at p < 0.05 out of 180 is fewer than the nine noise alone produces.

Those first three checks ask whether an association is internally consistent
given that you are looking at it. They do not ask whether you should have been
looking. The sweep-level correction asks that, and answers no.

greenyellow is the **extracellular matrix** module (interstitial matrix
q = 1.6e-06, ECM structural constituent q = 2.3e-05, ECM organization
q = 0.011; hubs LUM, FBLN2, ASPN, CALD1, TAGLN, MYH11). The direction is
mechanically coherent, since whole-muscle CSA includes interstitium while fibre
CSA does not, and the two are uncorrelated here at r = 0.03. Coherence is not
evidence, and this one did not clear the bar.

Its fit also rests on LR_S14, the highest subject on both axes — the same
subject V1's `06_mcsa_axis` result rested on, and a subject labelled Low
Responder despite the largest whole-muscle gain in the cohort.

### Biopsy composition is a real confound

At T2 the myofibre fraction correlates **−0.81** with change in whole-muscle
CSA and the blood fraction **+0.71**. Subjects who gained more muscle gave
less fibre-pure biopsies, which will induce proteome-wide differences tracking
that phenotype for reasons that are about the needle. Adjusting for the two
confounding panels removed four of the six within-cell hits, including
`HALLMARK_COAGULATION` — a blood signature that had appeared as a result.

## Reading the null

Nothing in this proteome tracks how much these subjects adapted, at any of
three feature levels, in any of six windows, against any of ten phenotypes. The
sweep returned fewer hits than chance predicts on every count, and no cell
survives correction across it.

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

`03_confirm_hits.R` permutes all 180 cells and takes about two hours.
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
