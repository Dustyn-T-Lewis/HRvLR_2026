# 02_Contrasts

Fits the nine contrasts with proteoDA and limma on the unimputed matrix.

| | |
|---|---|
| Reads | `01_Design/c_data/design.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/fit.rds`, `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf`, and proteoDA's per-contrast HTML reports in `b_reports/` (ignored) |
| Run | `Rscript 02_Differential_Expression/02_Contrasts/a_script/02_contrasts.R` |
| Cost | about 40 s |

`fit.rds` is the full proteoDA result, limma fit included; the pathway and module contrast steps
read it. BH runs within each contrast. `contrast_summary` gives tested, p < 0.05 and BH counts per
contrast. The fry concordance check reads the imputed matrix because fry takes no missing value.

The figure PDF opens with two overview pages. The biopsy overview shows all 1,900 proteins across
the 45 biopsies, centred per protein, with arm, time point, type I share, fCSA mixed, mCSA, leg
press 1RM and comp_hypertrophy above the columns. The biopsy traits are blank at T3, where none were
measured. The contrast overview shows each protein's log2 fold change in the nine contrasts beside
its status (not significant, p < 0.05, BH < 0.05, not tested). Rows follow average-linkage
clustering (1 - Pearson r) of the imputed matrix, for display only. After them come the p-value
histogram per contrast and a volcano per contrast, proteins at p < 0.05 in colour and proteins at BH
< 0.05 labelled. The log carries limma's note on partial NA coefficients for 34 proteins with an
empty cell.
