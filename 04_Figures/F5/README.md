# F5: modules

Co-expression modules from `01_Processing/06_Module_Building` (WGCNA, soft power 10, merge 0.25).
Colours are labels only, so every module is named by its top over-represented set (FDR < 0.05)
with the colour after it (`style/module_labels.R`). The per-module detail (hubs, rewiring,
coupling, per-person responses, PLS) is supplementary figure S7, below.

| Panel | Shows |
|---|---|
| a_heatmap | After YvO Figure 5A (`style/score_heatmap.R`): one row per module, largest first, with its size. Columns: the nine limma contrasts as r-equivalents of t (green–purple); then, for fibre area, muscle area and leg press, the T1 level against the T1 trait (baseline), the T1 level against the % gain (predicts) and the training change against the % gain (tracks), each as Spearman rho within HR and within LR (blue–red) and the HR − LR slope (green–purple). Outline: solid FDR < 0.05, dashed p < 0.05. |
| b_coupling | After the extra bout, the Spearman rho between every pair of modules' eigengene changes across the 8 people of each group (`02_HRvLR/08_Module_Detail`); boxed where HR and LR differ at Fisher p < 0.05. The training window is in `04_Figures/S7`. |
| c_preservation | `modulePreservation()` Zsummary for each group's own modules tested in the other group; dashed lines at 2 and 10, shading below 2. Each row names the group module and the full-cohort module holding most of its proteins. Point size is module size. |

## Inputs


| File | Used for |
|---|---|
| `01_Processing/06_Module_Building/c_data/module_building.rds` | eigengenes, biopsy table, module sizes |
| `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` | `labels` (top set), `string` |
| `02_HRvLR/03_Module_Contrasts/c_data/module_contrasts.rds` | `eigengene_tests` (a) |
| `02_HRvLR/08_Module_Detail/c_data/module_detail.rds` | `coupling` (b) |
| `02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds` | `within_group`, `group_difference` (a) |
| `02_HRvLR/05_Module_Preservation/c_data/module_preservation.rds` | `preservation`, `group_networks` |

## Read with care


- 6 of 90 eigengene tests reach p < 0.05 (4.5 expected) and none passes FDR < 0.05. The leads:
  aerobic respiration (red) and cristae formation (magenta) rise with training in HR and sit lower
  in HR than LR after the bout. S5 d draws the smallest p (magenta, HR − LR at T3, p 0.031) person by person.
- Three of 10 modules have no over-represented set at FDR < 0.05; they are named by colour alone.
- 7 of 10 modules re-form when any one person is left out (S5); green, yellow and magenta do not.
- Panel c shows most modules moderately preserved in the other group (Zsummary 2 to 10): HR and LR
  share their co-expression structure (26 of 29 group modules; 3 fall below 2, none above 10).
- In panel a, each within-group rho rests on 8 people: |rho| must reach about 0.71 for p < 0.05.
  The r-equivalent of a limma t uses its total df (410 for the eigengene contrasts), so a contrast at
  p < 0.05 shows near 0.1, and the HR − LR slope uses df 12.
- Panel b: after the bout, HR's modules change as one block and LR's do not. Respiration and
  cristae move together in HR (rho 0.90) and not in LR (0.17); mRNA processing, brown and
  cytoskeleton move together in HR and opposite to glycolysis. 6 of 45 pairs differ at Fisher
  p < 0.05. Each rho rests on 8 people.

## Supplementary figure S5: module characterisation and WGCNA QC

Supplement to F5. Modules are named as in F5 (`style/module_labels.R`).

| Panel | Shows |
|---|---|
| a_characterisation | Left: each module's three most over-represented sets (clusterProfiler ORA against the measured genes; GO Slim names that repeat a GO BP set counted once), x -log10 FDR capped at 20, fill the database, size the members in the set. Right: the module's ten highest-kME members, shaded by kME. |
| b_string | STRING edges (combined score ≥ 700) among each module's 25 highest-kME members, drawn with ggraph; each module laid out alone, components packed (`igraph::layout_components`). The five best-connected members are named. Strips give the module-wide STRING enrichment. |
| c_module_qc | WGCNA at soft powers 8, 12 and 14 against the power 10 used: modules found, proteins left in no module, and the share of each power-10 module's proteins kept together in its best match at the other power. |
| d_trajectories | Each person's eigengene at T1, T2 and T3, HR and LR, with the group means, for the module with the smallest eigengene p in a training, acute or HR − LR contrast (interactions left out); exploratory, the module is picked by p. |

### Inputs


| File | Used for |
|---|---|
| `01_Processing/06_Module_Building/c_data/module_building.rds` | module count and grey proteins at the power used |
| `01_Processing/06_Module_Building/c_data/module_building.xlsx` | `module_qc` |
| `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` | `labels`, `ora`, `top_members`, `string`, `network_edges` |

## Supplementary figure S7: module detail

Five PDFs from `02_HRvLR/08_Module_Detail`, drawn by `a_script/supp/S7_composite.R`. S7 draws the
step's tables whole, so its numbers are in `02_HRvLR/08_Module_Detail/c_data/module_detail.xlsx`
rather than a workbook here. Modules are named as in F5 (`style/module_labels.R`).

| File | Shows |
|---|---|
| `b_reports/supp/S7_hubs.pdf` | one page per module: its 15 highest-kME proteins (most central at the top) in the training and acute contrasts, HR, LR and HR − LR; fill logFC (capped at ±1), size −log10 p, filled at p < 0.05, black ring at FDR < 0.05 |
| `b_reports/supp/S7_rewiring.pdf` | one page per module: each member's kME within LR's 24 biopsies against within HR's; off the diagonal, a protein moves with the module in one group more than the other; the 8 largest shifts named |
| `b_reports/supp/S7_coupling.pdf` | Spearman rho between every pair of modules' eigengene changes (training, acute) within HR and within LR, and HR − LR boxed at Fisher p < 0.05 |
| `b_reports/supp/S7_per_person.pdf` | each person's eigengene change per module and window, people ordered by the growth score |
| `b_reports/supp/S7_pls.pdf` | sPLS-DA (HR against LR) and sPLS-R (trait gains) on the eigengenes: held-out score (balanced, repeated 4-fold CV) against its permutation null, and the component-1 loadings |

### Read with care

- Coupling: after the bout, HR's modules move as one block and LR's do not. Respiration (red)
  and cristae (magenta) change together in HR (rho 0.90) but not in LR (0.17); mRNA processing
  (blue), brown and cytoskeleton (green) change together in HR (0.88 to 0.93) and opposite to
  glycolysis (pink, -0.93 to -0.95), with no such pattern in LR. 7 of 90 pair-windows differ at
  Fisher p < 0.05, 6 of them after the bout. Each rho rests on 8 people.
- Rewiring: in the respiration module, NDUFA11 and the mitochondrially encoded MT-ND2 move with
  the module in HR (kME 0.73 and 0.82) and not in LR (-0.16 and 0.01); talin-2 (TLN2) does the same in
  the cristae module (0.78 against -0.21).
- PLS: no model beats its permutation null at p < 0.05. The closest is the training changes
  predicting the muscle-area gain (Q² 0.06, p 0.10); HR against LR is at chance from every window.
- No hub protein passes FDR < 0.05 in any contrast shown.

## Run

```sh
Rscript run_all.R 04_Figures/F5                       # F5 and S5, S7
Rscript 04_Figures/F5/a_script/F5_composite.R        # F5 alone
```

## Outputs

| File | Holds |
|---|---|
| `a_script/`, `b_reports/` | F5's panel scripts and PDFs (`a_*.R`, `a_*.pdf`, ...) and its composite (`F5_composite.pdf`: the figure, 250 × 185 mm landscape) |
| `a_script/supp/`, `b_reports/supp/` | S5: `S5_*.R`, `S5_*.pdf` and `S5_composite.pdf` (the figure, 250 × 210 mm: a beside c, then b beside d) |
| `a_script/supp/`, `b_reports/supp/` | S7: `S7_composite.R` and its five PDFs (`S7_hubs.pdf`, `S7_rewiring.pdf`, `S7_coupling.pdf`, `S7_per_person.pdf`, `S7_pls.pdf`) |
| `c_data/F5.xlsx`, `c_data/S5.xlsx` | one sheet per panel, the numbers drawn |
