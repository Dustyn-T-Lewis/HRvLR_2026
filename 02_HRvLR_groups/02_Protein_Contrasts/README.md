# 02_Protein_Contrasts

Which proteins change between time points, and which differ between HR and LR?

## The idea

A contrast is one comparison of group means. The design (`01_Design_Contrasts`) has six cell
means, one per group and time point (HR_T1, HR_T2, HR_T3, LR_T1, LR_T2, LR_T3), and nine contrasts
built from them:

| Family | Contrast | Compares |
|---|---|---|
| Training | `Training_HR`, `Training_LR` | T2 against T1 within a group |
| Acute | `Acute_HR`, `Acute_LR` | T3 against T2 within a group |
| HR - LR | `Baseline_HRvLR`, `Trained_HRvLR`, `Acute_HRvLR` | HR against LR at T1, T2, T3 |
| Interaction | `Training_Interaction`, `Acute_Interaction` | HR's change minus LR's change |

The interactions are the direct test of whether HR and LR respond differently. A protein that
changes in one group and not the other is not, on that basis, a group difference.
`Baseline_HRvLR` is the floor: HR and LR were defined from the outcome, so any difference at T1
shows how far the split itself separates people before training.

Each protein's logFC is the contrast on the log2 scale: +1 means twice as much in the later
time point (or in HR), -1 half as much. Its 95% confidence interval is the range of logFC the
data allow. The p-value asks how surprising the logFC would be if nothing changed; FDR
(Benjamini-Hochberg) adjusts p for the 1,900 proteins tested in that contrast, so FDR < 0.05 means
about 5% of the proteins called are expected to be false. Storey's q does the same with an
estimate of how many proteins have no effect.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Fit | `fit_limpa()` in `R/helpers.R`: `limpa::dpcDE()`, `limma::contrasts.fit()`, `limma::eBayes(robust = TRUE)` | limpa gives every protein a value in every biopsy and a standard error for each; uncertain values get less weight. Subject is a random block, because HR and LR are different people and each person gives up to three biopsies. Blood (haemoglobin index) is a covariate. |
| Effects | `limma::topTable(confint = TRUE)` | logFC with its 95% interval, moderated t, p and BH FDR per contrast; Storey's q from `qvalue`. |
| Blood check | the same fit without the blood term | Shows how much the blood covariate moves each result. |
| Null check | `limma::propTrueNull()`, median \|t\|, share of p below 0.2 | A calibrated null has median \|t\| near 0.67 and a flat p histogram. |
| Power | moderated standard error from `eBayes` | The logFC `Training_Interaction` could detect at 80% power, two-sided 0.05, per protein. |
| Concordance | `limma::fry()` | Takes one group's p < 0.05 proteins for a window and tests them as a set in the other group's contrast. The groups are different people, so this asks whether a response carries over without reusing the data that chose it. |
| Whole proteome | `prcomp()`, `vegan::adonis2()` | PCA on all proteins (scaled) and PERMANOVA on the distances between biopsies. Time is permuted within person; group is tested on one mean profile per person. |

## How to read the outputs

- Read logFC with its interval. With seven or eight people per group the intervals are wide; a
  nominal p < 0.05 with an interval crossing 0.1 says little about size.
- FDR is per contrast. A protein at FDR < 0.05 in one contrast and not another does not show the
  two contrasts differ.
- The concordance sets were chosen on p < 0.05, so their logFC in the choosing group is inflated;
  their logFC in the other group is not.
- PERMANOVA R² is the share of distance between biopsies the term explains; its p comes from
  permutations, not a model.

## Files

| | |
|---|---|
| Reads | `01_Design_Contrasts/c_data/design.rds` (design and contrasts), `01_Trunk/03_Quantification/c_data/proteins.rds` (limpa object), `01_Trunk/03_Quantification/c_data/DAList_corrected.rds` (blood-corrected complete matrix), `00_Input/phenotype.csv` |
| Writes | `c_data/protein_contrasts.xlsx`, `c_data/fit.rds` (the fit, read by `03_Set_Contrasts` and `04_Module_Contrasts`), `b_reports/protein_contrasts_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/02_Protein_Contrasts/a_script/protein_contrasts.R` |

```
02_Protein_Contrasts/
├── a_script/protein_contrasts.R
├── b_reports/protein_contrasts_report.pdf
└── c_data/
    ├── protein_contrasts.xlsx
    └── fit.rds
```

## Workbook: `protein_contrasts.xlsx`

`read_me` lists the sheets; `columns` defines every column.

| Sheet | One row is | Key columns |
|---|---|---|
| `protein_by_contrast` | a protein in a contrast | `logFC`, `logFC_lower`, `logFC_upper`, `t`, `p`, `fdr`, `qvalue`, `n_obs` |
| `contrast_summary` | a contrast | `n_tested`, `n_nominal`, `n_fdr_05`, `role` |
| `null_calibration` | a contrast | `prop_true_null`, `median_abs_t`, `p_below_0.2` |
| `detectable_effect` | a quantile | `log2_fc` detectable at 80% power |
| `fry_concordance` | a chosen set tested in another contrast | `set_from`, `ranked_on`, `p`, `fdr`, `mean_logfc_own`, `mean_logfc_other`, `same_sign_share` |
| `set_trajectory` | a set, group and time point | `mean_change` from each person's T1, `se` |
| `unadjusted` | a protein | logFC and p per contrast without the blood term |
| `pca`, `pca_variance` | a biopsy; a component | `PC1`, `PC2`; `variance_explained` |
| `permanova` | a test | `data`, `term`, `window`, `r2`, `f`, `p` |

## Report: `protein_contrasts_report.pdf`

1. Every protein across the 45 biopsies (centred), with group, time point and traits above.
2. logFC in the nine contrasts beside a p/FDR status map.
3. p-value histogram per contrast.
4. A volcano per contrast.
5. PCA of the biopsies.
