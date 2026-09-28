# 02_Contrasts

Fits the nine contrasts with proteoDA and limma on the unimputed matrix.

| | |
|---|---|
| Reads | `01_Design/c_data/design.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/fit.rds`, `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf`, and proteoDA's per-contrast HTML reports in `b_reports/` (ignored) |
| Run | `Rscript 02_Differential_Expression/02_Contrasts/a_script/02_contrasts.R` |
| Cost | about 40 s |

`fit.rds` is the full proteoDA result, limma fit included; the pathway and module contrast steps
read it. BH runs within each contrast. `contrast_summary` gives tested, p < 0.05 and BH counts per
contrast. The fry concordance check reads the imputed matrix because fry takes no missing value.

The figure PDF holds the p-value histogram per contrast and a volcano per contrast, proteins at
p < 0.05 in colour and proteins at BH < 0.05 labelled. The log carries limma's note on partial NA
coefficients for 34 proteins with an empty cell.
