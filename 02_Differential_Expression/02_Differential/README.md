# 02_Differential

Fits the model, writes the results and null checks, and classifies every protein on eight tasks.

| | |
|---|---|
| Reads | `01_Design/c_data/design.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/fit.rds`, `c_data/02_differential.xlsx`, `b_reports/02_differential_figures.pdf`, and proteoDA's nine `*_DA_report.html` with `static_plots/` in `b_reports/` (ignored) |
| Run | `quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports` |
| Cost | about 40 s |

`fit.rds` is the full proteoDA result, limma fit included. The figure PDF holds the p-value
histogram, then every protein nominal in any contrast and every protein nominal in any task as dot
matrices, 75 rows to a page.
