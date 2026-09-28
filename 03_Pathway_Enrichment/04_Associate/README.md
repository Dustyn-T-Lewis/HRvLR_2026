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

The figure PDF opens with four overview pages per collection: one for the sample model (within- and
between-person t for every trait, beside their status maps) and one per change-score window (pooled
rho for every outcome, beside its status map). Rows follow average-linkage clustering (1 - Pearson
r) of singscore within collection, for display only. The grids follow. Every model hit at p < 0.05
gets a panel: trait at the biopsy against singscore, one line per subject from T1 to T2, and a
header with the effect, p, BH q and the biopsies and subjects where the set was observed. Every
change-score pair at p < 0.05 gets a scatter, one point per subject, with a least-squares line per
arm drawn as a guide; the tests are Spearman. The header gives the pooled rho, p, BH q and subjects
used, the corner each arm's rho and exact p. Twelve panels to a page, ordered by p, then |t| or
|rho|, then ID. `plot_index` gives the PDF page and panel of every drawn row.
