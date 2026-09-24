# F01 · Cohort and filtering

| | |
|---|---|
| Reads | `00_Input/phenotype.csv`, `01_Preprocess/01_Filtering/c_data/01_filter.xlsx` |
| Writes | `b_reports/F01.pdf`, `b_reports/F01.png`, one PDF and PNG per panel in `b_reports/panels/`, `c_data/F01_data.xlsx` |
| Run | `Rscript 05_Figures/F01/a_script/F01.R` (each panel also runs alone) |

- A: four adaptation measures by arm, one point per subject.
- B: proteins left after each filter step.
- C: blood index by arm and timepoint, 45 analysed samples.

| File | Figure | Panel |
|---|---|---|
| `b_reports/F01.pdf`, `.png` | F01 | all |
| `b_reports/panels/A_phenotype.pdf`, `.png` | F01 | A |
| `b_reports/panels/B_filtering.pdf`, `.png` | F01 | B |
| `b_reports/panels/C_blood.pdf`, `.png` | F01 | C |
| `c_data/F01_data.xlsx`, sheet `A_phenotype` | F01 | A |
| `c_data/F01_data.xlsx`, sheet `B_filtering` | F01 | B |
| `c_data/F01_data.xlsx`, sheet `C_blood` | F01 | C |
