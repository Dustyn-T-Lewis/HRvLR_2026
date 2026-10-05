# F02_Groups

Figure 2: where HR and LR differ, using the group labels.

| | |
|---|---|
| Reads | `03_Pathways/03_Contrasts/c_data/set_tests.rds`, `03_Pathways/01_Gene_Sets/c_data/gene_sets.rds`, `03_Pathways/02_Scores/c_data/singscore.rds`, `04_Modules/01_Build/c_data/modules.rds`, `04_Modules/04_Associate/c_data/04_associate.xlsx` (`arm_difference`, `flips`), `05_Integration/03_Panel/c_data/03_panel.xlsx` (`panels`, `diablo`), the three `Classify` workbooks (`counts`), `00_Input/phenotype.csv` |
| Writes | `b_reports/f02_groups.pdf`, `b_reports/f02_groups.png` (300 x 110 mm, 600 dpi), `c_data/f02_groups.xlsx` |
| Run | `Rscript 06_Figures/F02_Groups/a_script/f02_groups.R` |

The script reads saved outputs only and checks its plotted numbers against them before writing.

## Panels

| Panel | Shows | Reading |
|---|---|---|
| a | The six strongest fgsea sets per interaction contrast (HR change minus LR change), overlapping sets dropped | mitochondrial sets rise more in HR over training and fall more in HR after the bout; Reactome aerobic respiration NES -2.98 acute |
| b | Each person's Hallmark OXPHOS score at T1, T2, T3, arm means bold | HR starts lower (HR - LR NES -2.60 at T1), reaches LR at T2 (training interaction NES +2.25) and drops after the bout |
| c | Each person's T3 - T2 change in that score | HR falls (median -0.012), LR holds (median about 0); interaction NES -2.70, but fry p = 0.15 |
| d | Black (sarcomere) module training change against leg-press gain, one fitted line per arm | slopes run in opposite directions: rho HR -1.00, LR +0.48; arm x gain FDR 0.0045 |
| e | Cross-validated HR vs LR AUC for the module and protein panels and the DIABLO models, against 200 label permutations | every model sits at or below chance; no protein or pathway separates the arms at FDR < 0.05 (modules: 2) |

The OXPHOS results rest on fgsea and camera, which rank proteins; fry, which rotates samples and so
keeps the correlation between proteins, does not support them. Panel d has seven or eight people per
arm. The arms differ in how they change, not in a separable signature.

## Workbook

| Sheet | Holds |
|---|---|
| `interaction_sets` | Panel a: plotted sets per interaction contrast, with NES and fgsea FDR. |
| `oxphos_scores` | Panel b: each sample's Hallmark OXPHOS singscore. |
| `acute_change` | Panel c: each person's T3 - T2 change. |
| `black_change` | Panel d: each person's black module training change and leg-press gain. |
| `classification` | Panel e: AUC and permutation null per panel and DIABLO model. |
| `single_features` | Single-feature HR vs LR tests at FDR < 0.05, per level. |
