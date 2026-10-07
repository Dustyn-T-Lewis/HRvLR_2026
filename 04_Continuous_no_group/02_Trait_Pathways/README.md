# 02_Trait_Pathways

Which gene sets lead each trait's protein ranking, every subject pooled.

| | |
|---|---|
| Reads | `01_Trait_Associations/01_Proteins/c_data/associate.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/02_pathways.xlsx`, `c_data/pathways.rds`, `b_reports/02_pathways_figures.pdf` |
| Run | `Rscript 04_Continuous_no_group/02_Trait_Pathways/a_script/02_pathways.R` |

## Method

Per pairing and trait, proteins are ranked by the limma t from `01_Trait_Associations/01_Proteins`. fgsea
(seeded) tests every set on that ranking; its p and BH are its own.

fgsea permutes proteins, and proteins in one set move together, so set counts here run high on
any ranking. Read the leading sets, not the counts.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Per window and trait: sets tested, at p < 0.05, at FDR < 0.05 (up, down). |
| `pathways` | Per window, trait and set: fgsea NES, p and BH. |
