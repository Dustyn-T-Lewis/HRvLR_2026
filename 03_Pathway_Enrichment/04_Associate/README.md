# 04_Associate

Ties each set's singscore to phenotype, the protein level's two ways.

| | |
|---|---|
| Reads | `01_Scores/c_data/singscore.rds`, `00_Gene_Sets/c_data/gene_sets.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` (sample sheet), `00_Input/phenotype.csv` |
| Writes | `c_data/04_associate.xlsx`, `b_reports/04_associate_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/04_Associate/a_script/04_associate.R` |
| Cost | about 3 min |

The sample model and the change score are those of `02_Differential_Expression/04_Associate`, on
singscore instead of abundance. Scores have no missing values, so no observation floor applies. BH
runs within trait, term and collection (model) or window, outcome and collection (change score).
`counts` gives tested, p < 0.05 and BH < 0.05 counts per collection.

The figure PDF gives every model hit at p < 0.05 a panel (one line per subject from T1 to T2),
then every change-score pair at p < 0.05 a scatter (`ggpubr::ggscatter`, a line per arm, each
arm's rho and exact p in the corner), twelve to a page, by p.
