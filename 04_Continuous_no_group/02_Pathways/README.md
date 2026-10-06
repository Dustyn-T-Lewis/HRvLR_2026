# 02_Pathways

Which gene sets lead each trait's protein ranking, every subject pooled.

| | |
|---|---|
| Reads | `01_Associate/01_Proteins/c_data/associate.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/02_pathways.xlsx`, `c_data/pathways.rds`, `b_reports/02_pathways_figures.pdf` |
| Run | `Rscript 04_Continuous_no_group/02_Pathways/a_script/02_pathways.R` |

## Method

Per pairing and trait, proteins are ranked by the limma t from `01_Associate/01_Proteins`. fgsea
(seeded) tests every set on that ranking; its p and BH are its own. cameraPR runs beside it, and
`camera_agrees` marks sets fgsea calls at FDR < 0.05 that cameraPR reaches at p < 0.05 in the same
direction. cameraPR uses limma's preset inter-gene correlation (0.01), so it is close to fgsea, not
an independent check.

fgsea permutes proteins, and proteins in one set move together, so set counts here run high on
any ranking. Read the leading sets, not the counts.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Per window and trait: sets tested, at p < 0.05, at FDR < 0.05 (up, down), cameraPR agreeing. |
| `pathways` | Per window, trait and set: NES, p, BH; cameraPR direction and p; `camera_agrees`. |
