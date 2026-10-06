# 01_Build

Builds WGCNA modules and scores their eigengenes.

| | |
|---|---|
| Reads | `01_Trunk/03_Quantify/c_data/DAList_corrected.rds` |
| Writes | `c_data/modules.rds`, `c_data/01_build.xlsx`, `b_reports/01_build_figures.pdf` |
| Run | `Rscript 01_Trunk/06_Modules/a_script/06_modules.R` |

## Method

WGCNA on subject-centred abundance (limpa's complete matrix with the blood slope removed), so modules do not encode who a biopsy came from;
eigengenes are scored on raw abundance. A subject with one biopsy centres to zeros and is left out
of the definition. Signed network, biweight midcorrelation (`maxPOutliers = 0.05`), `deepSplit = 2`,
`minModuleSize = 30`, `mergeCutHeight = 0.15`; the soft power is `pickSoftThreshold()`'s estimate at
R2 0.85. Membership (kME) is read on the same subject-centred matrix that defined the modules. The
subject ICC fits `eigengene ~ timepoint + (1 | subject)`, so time effects do not count as noise.
`modules.rds` feeds steps 02 to 06.

## Figures

Soft-threshold curves, module sizes with subject ICC, and eigengene trajectories by arm.

## Workbook

| Sheet | Holds |
|---|---|
| `module_summary` | Proteins per module and the subject ICC of its eigengene, time point fitted. |
| `membership` | Every protein's module and its kME in that module, on the subject-centred matrix. |
| `eigengenes` | Module eigengene per sample, scored on raw abundance. |
| `soft_threshold` | pickSoftThreshold fit indices per power; chosen marks the power used. |
