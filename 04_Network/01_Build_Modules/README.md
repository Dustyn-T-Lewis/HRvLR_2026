# 01_Build_Modules

Builds WGCNA modules and scores their eigengenes.

| | |
|---|---|
| Reads | `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/modules.rds`, `c_data/01_build_modules.xlsx`, `b_reports/01_build_modules_figures.pdf` |
| Run | `Rscript 04_Network/01_Build_Modules/a_script/01_build_modules.R` |

## Method

WGCNA on subject-centred imputed abundance, so modules do not encode who a biopsy came from;
eigengenes are scored on raw abundance. A subject with one biopsy centres to zeros and is left out
of the definition. Signed network, biweight midcorrelation (`maxPOutliers = 0.05`), `deepSplit = 2`,
`minModuleSize = 30`, `mergeCutHeight = 0.15`; the soft power is `pickSoftThreshold()`'s estimate at
R2 0.85. `modules.rds` feeds steps 02 to 06.

## Figures

Soft-threshold curves, module sizes with subject ICC, and eigengene trajectories by arm.

## Workbook

| Sheet | Holds |
|---|---|
| `module_summary` | Proteins per module and the subject ICC of its eigengene. |
| `membership` | Every protein's module and its kME in that module. |
| `eigengenes` | Module eigengene per sample, scored on raw abundance. |
| `soft_threshold` | pickSoftThreshold fit indices per power; chosen marks the power used. |
