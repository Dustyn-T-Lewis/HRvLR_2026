# 02_Differential_Expression / 03_Phenotype

Correlates each protein with the ten adaptation measures.

| | |
|---|---|
| Script | `a_script/03_phenotype.qmd` |
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/phenotype.rds`, `c_data/03_phenotype.xlsx`, `b_reports/03_phenotype_figures.pdf` |

A contrast compares group means and cannot use a per-subject outcome, so this is a separate
question. Each subject contributes one value per protein in three windows: training change
(T2 − T1, the window the phenotypes share), baseline level (T1) and acute change
(T3 − T2). Spearman, t approximation; BH within window and outcome. A protein is tested only when
two thirds of the window's subjects observed it.

| Window | Outcome | Nominal / chance | BH < 0.05 |
|---|---|---:|---:|
| training | d_1rm_ext | 1.96 | 0 |
| training | d_mcsa | 0.95 | 1 (RPS4X, rho 0.91, n 14) |
| baseline | d_1rm_ext | 2.36 | 0 |
| baseline | d_1rm_legpress | 1.25 | 1 (ACSL3, rho −0.87, n 15) |
| baseline | volume_load | 1.46 | 0 |
| acute | comp_hypertrophy | 1.86 | 0 |
| acute | d_mcsa | 1.52 | 0 |

The other 23 cells sit between 0.44 and 1.18. The full table is the workbook's `summary` sheet.
With 14 or 15 subjects a correlation needs to reach about 0.7 to be detected at 80% power, so the
screen is null at that resolution apart from RPS4X and ACSL3.

`cor.test`'s "exact" Spearman p is an Edgeworth series above nine observations and returns 0 in
the far tail, which put two proteins at FDR 0. The t approximation replaces it.

Figures: nominal hits over chance by outcome and window, then one hit matrix per window.
`by_arm` in the workbook repeats the correlations inside each arm and window, descriptively; each
correlation rests on 4 to 8 subjects.
