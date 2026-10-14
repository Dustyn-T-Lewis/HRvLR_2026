# 06_Module_Building

Finds co-expression modules with WGCNA on person-centred abundance, then scores each module's
eigengene on the uncentred matrix so differences between people stay testable. No label or
phenotype enters construction.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/03_Quantification/c_data/DAList_corrected.rds` | limpa's complete matrix with blood slopes removed |

## Run

```sh
quarto render 01_Processing/06_Module_Building/a_script/module_building.qmd --output-dir ../b_reports
```

Run from the project root, so R starts with the project's renv library.

## Outputs

| File | Holds |
|---|---|
| `b_reports/module_building.html` | the walkthrough: centring, soft power, modules, kME, person ICC, trajectories |
| `c_data/modules.rds` | eigengenes, membership with kME, module summary; read by every module-level step |
| `c_data/module_building.xlsx` | the same tables and the soft-power fit indices |

Module colours are labels only. They change if `01_Processing` is rebuilt.
