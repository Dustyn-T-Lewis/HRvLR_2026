# 01_Permutation

Gives the candidate counts a chance line: how many features reach k nominal hits when the link
between a person's proteome and their response is broken. Context only; nothing is filtered.

| | |
|---|---|
| Reads | `01_Preprocess/03_Quantification/c_data/proteins.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds`, `03_Pathway_Enrichment/01_Scores/c_data/singscore.rds`, `04_Network/01_Build_Modules/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/01_permutation.xlsx`, `b_reports/01_permutation_figures.pdf` |
| Run | `Rscript 05_Summary/01_Permutation/a_script/01_permutation.R` (about 15 minutes on eight cores) |

## Method

Each of 200 shuffles hands every subject another subject's phenotype row, arm and traits together,
so arm sizes and the correlations between traits hold. The tests behind the candidate lists are
rerun on the shuffled labels at all three levels: the HR-against-LR contrasts (limpa for proteins,
fgsea and fry for sets, eigengene limma and fry for modules), ANCOVA, the HR-against-LR rank tests,
and the training, baseline and acute associations. Blood correction is redone with the shuffled
groups. Hits are counted with `on_lists()` in `R/helpers.R`, the rules `02_Candidates` uses.

The real labelling runs through the same code as shuffle 0, so observed and permuted counts share
one code path. Two simplifications keep 200 shuffles tractable: within-arm Spearman p uses the t
approximation rather than the exact permutation p (so this step's observed counts run slightly
above `02_Candidates`'), and set scores and module eigengenes are taken as computed once.

Per feature, `permutation_score` is the share of shuffles in which it reached at least its observed
hit count, (b + 1)/(m + 1) (Phipson & Smyth 2010). It is more lenient than a per-test correction
because it keeps the correlation between tests, and it is exploratory context, not a gate.

## Figures

One page: per level and list, features with at least k hits under the real labelling against the
permuted median and 95th percentile.

## Workbook

| Sheet | Holds |
|---|---|
| `chance_line` | Per level, list and k: real count, permuted median and 95th percentile, share of shuffles at or above the real count. |
| `feature_scores` | Per feature: hits under the real labelling and the share of shuffles reaching that many. |
