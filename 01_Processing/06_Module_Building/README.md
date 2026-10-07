# 06_Module_Building

Which proteins rise and fall together within a person, and how does each such group move in each biopsy?

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/DAList_corrected.rds` |
| Writes | `c_data/module_building.xlsx`, `c_data/modules.rds`, `b_reports/module_building_report.pdf` |
| Run | `Rscript 01_Processing/06_Module_Building/a_script/module_building.R` |

## Key result

- 15 modules of 41 to 285 proteins; 304 proteins fit no module (`grey`).
- Soft power 11 (scale-free fit R² 0.86).

## Method

- WGCNA (`blockwiseModules()`): signed network, bicor, `deepSplit = 2`, `minModuleSize = 30`, `mergeCutHeight = 0.15`, `randomSeed = 54321`; soft power from `pickSoftThreshold()` (R² 0.85).
- Modules are defined on abundance centred within each person, so they describe co-movement within people rather than differences between them.
- Eigengenes (`moduleEigengenes()`) are scored on uncentred abundance in all 45 biopsies; kME (`signedKME()`) is each protein's correlation with its module's eigengene. `grey` holds proteins in no module.
- `icc` is how much of an eigengene's variance lies between people (`lme4::lmer`).

## Read with care

- Module colours are names, not ranks, and change when 01_Processing is rebuilt; compare modules by membership.

## Sheets

`module_building.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `module_summary` | Proteins per module and the subject ICC of its eigengene. |
| `membership` | Every protein's module and its kME in that module. |
| `eigengenes` | Module eigengene per sample, scored on raw abundance. |
| `soft_threshold` | pickSoftThreshold fit indices per power; chosen marks the power used. |
