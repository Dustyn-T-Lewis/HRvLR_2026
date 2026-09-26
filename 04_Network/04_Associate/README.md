# 04_Associate

Ties each module eigengene to phenotype, the protein level's two ways.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/04_associate.xlsx`, `b_reports/04_associate_figures.pdf` |
| Run | `Rscript 04_Network/04_Associate/a_script/04_associate.R` |
| Cost | about 10 s |

The sample model and the change score are those of `02_Differential_Expression/04_Associate`, on
eigengenes instead of abundance. Eigengenes have no missing values, so no observation floor
applies. BH runs within trait and term (model) or window and outcome (change score), over twelve
modules. `change_score_by_arm` keeps every within-arm correlation with its exact p.

The figure PDF opens on the chance ratios, then gives every nominal model hit a panel (one line per
subject from T1 to T2), then every nominal change-score pair a scatter (`ggpubr::ggscatter`, a line
per arm, each arm's rho and exact p in the corner).
