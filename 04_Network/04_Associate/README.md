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

The figure PDF opens with overview pages: one for the sample model (within- and between-person t for
every trait, beside their status maps) and one per change-score window (pooled rho for every
outcome, beside its status map). Rows follow module size. The grids follow. Every model hit at p <
0.05 gets a panel: trait at the biopsy against the eigengene, one line per subject from T1 to T2,
and a header with the effect, p, BH q and the biopsies and subjects where the module was observed.
Every change-score pair at p < 0.05 gets a scatter, one point per subject, with a least-squares line
per arm drawn as a guide; the tests are Spearman. The header gives the pooled rho, p, BH q and
subjects used, the corner each arm's rho and exact p. Twelve panels to a page, ordered by p, then
|t| or |rho|, then ID. `plot_index` gives the PDF page and panel of every drawn row.
