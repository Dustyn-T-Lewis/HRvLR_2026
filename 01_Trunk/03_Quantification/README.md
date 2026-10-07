# 03_Quantification

How do we give every protein a value in every sample without pretending the missing ones were
measured?

## The idea

In data-independent acquisition (DIA) proteomics a protein is missing from a run mostly because it
was too faint to detect there. A missing value therefore carries information: the protein was
probably low. Dropping it throws that away; filling it with a constant or a guess treats an
estimate as a measurement and makes the protein look more certain than it is.

limpa (Li, Cobbold & Smyth 2025) models the missingness instead. It fits one detection probability
curve (DPC) for the whole experiment: the probability that a value is detected, as a logistic
function of its true log2 intensity. Two numbers describe it, an intercept and a slope; together
they place the intensity at which half the values are detected, and the slope sets how sharply
detection rises with intensity. A steep slope means missing values sit just below the detection
limit; a shallow one means detection is noisier and a missing value says less.

With the curve in hand, limpa gives each protein a value in every sample (`dpcQuantByRow`). An
observed value stays close to what was measured, with a small standard error. A missing value gets
the value expected from the curve and the protein's observed values, with a wide standard error.
The protein-level fits downstream (`limpa::dpcDE`) turn the standard errors into precision weights,
so a protein estimated in many samples counts for less than one measured in all of them.

Our input is already summarised to proteins, so each row is its own feature: the route the limpa
vignette gives for protein-level input. With precursor-level input limpa would also pool peptides;
here there are none to pool.

Blood is handled in two ways, because not every method can carry a covariate. limma models carry
the blood index as a covariate. Methods that cannot (rank scores, WGCNA, correlations) read a copy
with each protein's blood slope removed. The slope is fitted beside time point only
(`limma::removeBatchEffect`, design `~ 0 + timepoint`), so no group label shapes the corrected
matrix. The index is centred to mean zero first, so the correction moves no protein's average.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Curve | `limpa::dpc()` | One DPC for the run, from the normalised log2 matrix. |
| Quantify | `limpa::dpcQuantByRow()` | A complete matrix `E` with a standard error for every value; the script checks that no value is missing and that the observed mask matches the input. |
| Blood | `blood_corrected()` in `R/helpers.R`: `limma::removeBatchEffect(covariates = blood, design = ~ 0 + timepoint)` | `blood` is `blood_rel` from `01_Protein_Filtering`, centred. |
| Sample correlation | `cor()` | Pearson correlation between samples on the corrected matrix, each protein centred on its mean, so it reflects how samples differ rather than the abundance order every sample shares. |
| Report | `limpa::plotDPC()` | The fitted curve with its intercept, slope and share of values missing. |

## How to read the outputs

- `curve` gives the DPC on the logit scale. A higher slope means missingness tracks intensity more
  tightly, so each missing value is more informative.
- In `proteins`, `prop_observed` says how much of a protein's row was measured. A protein with a
  low share rests on the curve more than on data; its `median_se_missing` shows how wide those
  estimates are.
- `mean_observed` averages the measured values only; `mean_complete` averages the whole row,
  estimates included.
- `sample_correlation` uses the blood-corrected matrix. `01_Trunk/Figures/S01_QC` summarises it
  within and between people; each person gives up to three biopsies, which is why the arcs block
  on subject.

## Files

| | |
|---|---|
| Reads | `01_Trunk/02_Normalisation/c_data/DAList_normalized.rds` |
| Writes | `c_data/proteins.rds`, `c_data/DAList_limpa.rds`, `c_data/DAList_corrected.rds`, `c_data/quantification.xlsx`, `b_reports/quantification_report.pdf` |
| Run | `Rscript 01_Trunk/03_Quantification/a_script/quantification.R` |

```
03_Quantification/
├── a_script/quantification.R
├── b_reports/quantification_report.pdf
└── c_data/
    ├── quantification.xlsx
    ├── proteins.rds
    ├── DAList_limpa.rds
    └── DAList_corrected.rds
```

The three `.rds` files:

| File | Holds | Read by |
|---|---|---|
| `proteins.rds` | `proteins` (limpa's EList: complete values, standard errors, observation counts), `metadata` (sample sheet with `group` and centred `blood`), `annotation` | protein contrasts (`fit_limpa()` in `R/helpers.R`) and the HR/LR design; copied and relabelled by each `03_Cluster_groups/*/01_Clustering` |
| `DAList_limpa.rds` | `data` (the complete matrix, no blood correction), `annotation`, `metadata`, `observed` (TRUE where the value was measured) | fry in the set and module contrasts (blood in the design), row order of the overview heatmaps in protein classification and protein-trait associations; copied by each `01_Clustering` |
| `DAList_corrected.rds` | the same, with each protein's blood slope removed | `05_Set_Scores`, `06_Module_Building`, and the rank tests, correlations, PCA and figures in the arcs; copied by each `01_Clustering` |

`DAList_limpa.rds` and `DAList_corrected.rds` are plain lists shaped like a proteoDA DAList, with
the extra `observed` mask.

## Workbook: `quantification.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close it.

| Sheet | One row is | Key columns |
|---|---|---|
| `curve` | a DPC term (intercept, slope) | `estimate`, `pct_missing` |
| `proteins` | a protein | `n_observed`, `prop_observed`, `mean_observed`, `mean_complete`, `median_se_observed`, `median_se_missing` |
| `blood` | a sample | `subject`, `group`, `blood_rel`, `blood` (centred) |
| `sample_correlation` | a sample | Pearson r with every other sample |

## Report: `quantification_report.pdf`

One page: the detection probability curve, with its intercept, slope and the share of values
missing in the subtitle.
