# F5: modules

Co-expression modules from `01_Processing/06_Module_Building` (WGCNA, soft power 10, merge 0.25).
Colours are labels only, so every module is named by its top over-represented set (FDR < 0.05)
with the colour after it (`style/module_labels.R`). The per-module detail (hubs, rewiring,
coupling, per-person responses, PLS) is in `04_Figures/Modules`.

| Panel | Shows |
|---|---|
| a_heatmap | After YvO Figure 5A (`style/score_heatmap.R`): one row per module, largest first, with its size. Columns: the nine limma contrasts as r-equivalents of t (green–purple); then, for fibre area, muscle area and leg press, the T1 level against the T1 trait (baseline), the T1 level against the % gain (predicts) and the training change against the % gain (tracks), each as Spearman rho within HR and within LR (blue–red) and the HR − LR slope (green–purple). Outline: solid FDR < 0.05, dashed p < 0.05. |
| b_trajectories | Each person's eigengene at T1, T2 and T3, HR and LR, with the group means, for the module with the smallest eigengene p in a training, acute or HR − LR contrast (interactions left out). |
| c_preservation | `modulePreservation()` Zsummary for each group's own modules tested in the other group; dashed lines at 2 and 10, shading below 2. Each row names the group module and the full-cohort module holding most of its proteins. Point size is module size. |

## Inputs

| File | Used for |
|---|---|
| `01_Processing/06_Module_Building/c_data/modules.rds` | eigengenes, biopsy table, module sizes |
| `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` | `labels` (top set), `string` |
| `02_HRvLR/03_Module_Contrasts/c_data/module_tests.rds` | `eigengene_tests` (a, b) |
| `02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds` | `within_group`, `group_difference` (a) |
| `02_HRvLR/05_Module_Preservation/c_data/preservation.rds` | `preservation`, `group_networks` |

## Run

```sh
Rscript 04_Figures/F5/a_script/F5_composite.R                 # figure and workbook
Rscript 04_Figures/F5/a_script/panels/a_heatmap.R             # one panel alone
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/panels/<panel>.pdf` | each panel |
| `b_reports/F5_composite.pdf` | the figure, 250 × 180 mm landscape |
| `c_data/F5.xlsx` | one sheet per panel with the numbers drawn |

## Read with care

- 6 of 90 eigengene tests reach p < 0.05 (4.5 expected) and none passes FDR < 0.05. The leads:
  aerobic respiration (red) and cristae formation (magenta) rise with training in HR and sit lower
  in HR than LR after the bout. Panel b draws the smallest p (magenta, HR − LR at T3, p 0.031).
- Three of 10 modules have no over-represented set at FDR < 0.05; they are named by colour alone.
- 7 of 10 modules re-form when any one person is left out (S5); yellow and magenta do not.
- Panel c shows most modules moderately preserved in the other group (Zsummary 2 to 10): HR and LR
  share their co-expression structure (26 of 29 group modules; 3 fall below 2, none above 10).
- In panel a, each within-group rho rests on 8 people: |rho| must reach about 0.71 for p < 0.05.
  The r-equivalent of a limma t uses its total df (about 46 for contrasts), so a contrast at
  p < 0.05 shows near 0.3, and the HR − LR slope uses df 12.
