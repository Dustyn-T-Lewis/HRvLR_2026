# 01_Build_Modules

Builds WGCNA modules and scores their eigengenes.

| | |
|---|---|
| Reads | `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/modules.rds`, `c_data/01_build_modules.xlsx`, `b_reports/01_build_modules_figures.pdf` |
| Run | `Rscript 04_Network/01_Build_Modules/a_script/01_build_modules.R` |
| Cost | about 10 s |

`modules.rds` holds `eigengenes`, `membership` (module and kME per protein), `meta` and
`module_summary`, which steps 02 to 06 read. The figure PDF holds the soft-threshold curves, module
sizes with subject ICC, and eigengene trajectories by arm.
