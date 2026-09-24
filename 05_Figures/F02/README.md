# F02 · Proteome

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/02_normalize.xlsx`, `02_Differential_Expression/02_Differential/c_data/02_differential.xlsx` |
| Writes | `b_reports/F02.pdf`, `b_reports/F02.png`, one PDF and PNG per panel in `b_reports/panels/`, `c_data/F02_data.xlsx` |
| Run | `Rscript 05_Figures/F02/a_script/F02.R` (each panel also runs alone) |

- A: samples on PC1 and PC2 after normalisation.
- B: nominal proteins per contrast against 5% of those tested.
- C: p-value histograms for Training_Interaction and Baseline_HRvLR.
- D: Training_Interaction, every tested protein.

| File | Figure | Panel |
|---|---|---|
| `b_reports/F02.pdf`, `.png` | F02 | all |
| `b_reports/panels/A_pca.pdf`, `.png` | F02 | A |
| `b_reports/panels/B_hits.pdf`, `.png` | F02 | B |
| `b_reports/panels/C_p_histogram.pdf`, `.png` | F02 | C |
| `b_reports/panels/D_volcano.pdf`, `.png` | F02 | D |
| `c_data/F02_data.xlsx`, sheet `A_pca` | F02 | A |
| `c_data/F02_data.xlsx`, sheet `B_hits` | F02 | B |
| `c_data/F02_data.xlsx`, sheet `C_p_histogram` | F02 | C |
| `c_data/F02_data.xlsx`, sheet `D_volcano` | F02 | D |
