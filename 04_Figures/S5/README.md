# S5: module characterisation and WGCNA QC

Supplement to F5. Modules are named as in F5 (`style/module_labels.R`).

| Panel | Shows |
|---|---|
| a_characterisation | Left: each module's three most over-represented sets (clusterProfiler ORA against the measured genes; GO Slim names that repeat a GO BP set counted once), x -log10 FDR capped at 20, fill the database, size the members in the set. Right: the module's ten highest-kME members, shaded by kME. |
| b_string | STRING edges (combined score ≥ 700) among each module's 25 highest-kME members, drawn with ggraph; each module laid out alone, components packed (`igraph::layout_components`). The five best-connected members are named. Strips give the module-wide STRING enrichment. |
| c_module_qc | WGCNA at soft powers 10, 14 and 19 against the power 12 used: modules found, proteins left in no module, and the share of each power-12 module's proteins kept together in its best match at the other power. |

## Inputs

| File | Used for |
|---|---|
| `01_Processing/06_Module_Building/c_data/modules.rds` | power-12 module count and grey proteins |
| `01_Processing/06_Module_Building/c_data/module_building.xlsx` | `module_qc` |
| `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` | `labels`, `ora`, `top_members`, `string`, `network_edges` |

## Run

```sh
Rscript 04_Figures/S5/a_script/S5_composite.R              # figure and workbook
Rscript 04_Figures/S5/a_script/panels/b_string.R           # one panel alone
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel |
| `b_reports/S5_composite.pdf` | the figure, 180 × 238 mm |
| `c_data/S5.xlsx` | one sheet per panel with the numbers drawn |
