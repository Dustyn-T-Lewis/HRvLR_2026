# 03_Phenotype

Correlates each protein with the ten phenotypes in the training, baseline and acute windows.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/phenotype.rds`, `c_data/03_phenotype.xlsx`, `b_reports/03_phenotype_figures.pdf` |
| Run | `quarto render 02_Differential_Expression/03_Phenotype/a_script/03_phenotype.qmd --output-dir ../b_reports` |
| Cost | about 25 s |

`phenotype.rds` holds `protein_association`, which `04_Network/04` reads. The figure PDF holds
nominal hits over chance by outcome and window, then one hit matrix per window.
