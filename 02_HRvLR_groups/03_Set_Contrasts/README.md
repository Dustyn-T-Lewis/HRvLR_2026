# 03_Set_Contrasts

Which gene sets shift in each of the nine contrasts, judged two ways?

## The idea

A gene set is a published list of genes that share a function or pathway: oxidative
phosphorylation, the ribosome, muscle contraction. The sets come from five collections (Hallmark,
Reactome, KEGG_Legacy, GOBP, GO_Slim) and are mapped to proteins in `01_Trunk/04_Gene_Sets`, one
representative protein per gene symbol. A set is tested when 15 to 500 of its proteins were
measured. Testing a set rather than a protein asks whether a function moved, and pools many small
protein changes into one question.

Two methods run on every set and contrast. They ask different questions, so they can disagree.

fgsea is competitive: does this set rank above or below the other proteins? It sorts every protein
by its moderated t in the contrast (from `02_Protein_Contrasts`), then walks down the list, adding
to a running score at each member and subtracting at each non-member. The enrichment score is the
largest distance from zero that walk reaches. The normalised enrichment score (NES) divides it by
the average score of random sets of the same size, so sets of different sizes compare: positive
means the set sits toward the top (up in the contrast), negative toward the bottom. The p-value
asks how often a random set of proteins of the same size scores as far out. That null shuffles
proteins, not people, so it treats each protein as an independent draw. Proteins in a pathway are
correlated, so a large, tightly co-regulated set can reach a small p with few people behind it.
The leading edge is the members the walk passes before it reaches that peak: the proteins driving
the score.

fry is self-contained: did this set's proteins move away from zero, together, in these people? It
works inside the linear model (same six cell means, blood slope, subject block), averages the
members' evidence for the contrast, and compares it with rotations of the residuals. Rotation
keeps the correlation between proteins, so 150 OXPHOS proteins rising in the same few people count
as one consistent signal, not 150 confirmations. Its p therefore reflects how many people carry
the shift, and with seven or eight per group that bar is high. fry reports a direction, Up or
Down; the stage keeps its directional p.

They disagree for known reasons. If most of the proteome moves up, a set that moves with it is not
enriched against the rest (fgsea) but has still moved (fry). A correlated set gets a small fgsea p
from its size alone; fry does not count its members as independent.

Sets overlap: glycolysis appears in Hallmark, KEGG, Reactome and GO. fgsea's `collapsePathways()`
takes the significant sets strongest first and re-tests each one conditioned on a stronger set's
leading edge; a set that stays significant is kept as non-redundant (`main = TRUE`). It marks
sets and deletes none.

FDR (Benjamini-Hochberg) is the p-value adjusted for every set the method tested in that contrast,
all five collections pooled; at FDR < 0.05 about 5% of the sets called are expected to be false.
Storey's q does the same with an estimate of how many sets have no effect.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Protein results | `limma::topTable()` on the saved fit | logFC, t, p and BH FDR per protein and contrast, joined to the gene map. `pi_score` = p^\|logFC\| orders proteins for figures and selects nothing. |
| Ranking | one named vector per contrast | Representative proteins only, keyed by gene symbol so leading edges come back as symbols. |
| fgsea | `fgsea::fgsea(minSize = 15, maxSize = 500)` | Seeded (`set.seed(1)`); NES, p and fgsea's own BH per contrast. |
| fry | `limma::fry(block = subject, correlation = ...)` | On limpa's complete matrix (fry takes no missing value) with the design from `design.rds`. The within-subject correlation is re-estimated on that matrix with `limma::duplicateCorrelation()`. Gene sets are mapped back to matrix rows through the representative protein. |
| Storey q | `storey_q()` in `R/helpers.R` | Within contrast and method. |
| Redundancy | `fgsea::collapsePathways(pval.threshold = 0.05)` | On the fgsea sets at FDR < 0.05; sets the `main` flag. |
| Overview scores | singscore from `01_Trunk/05_Set_Scores` | For the biopsy heatmap only. singscore places a set within each biopsy's own ranking of proteins, so every biopsy gets one score per set. |

## How to read the outputs

- Read the floor first. `set_summary` lists `Baseline_HRvLR` at the top: the calls the HR/LR split
  produces before training. A contrast's count means little until set against it.
- An fgsea call without a fry call says the set ranks high against other proteins, not that the
  people behind it moved consistently. A fry call says the latter.
- fry needs every member present, so it runs on limpa's filled-in matrix; fgsea ranks on the
  limpa fit, which weights each value by its precision. The two do not see identical inputs.
- Each method's FDR is per contrast. A set at FDR < 0.05 in `Training_HR` and not in
  `Training_LR` does not show the groups differ; the interaction contrasts test that.
- `collapsePathways()` prunes after BH, so the FDR of a kept set is conservative, and a pruned set
  is not a false one, only overlapping.

## Files

| | |
|---|---|
| Reads | `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `02_HRvLR_groups/02_Protein_Contrasts/c_data/fit.rds`, `02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds`, `01_Trunk/03_Quantification/c_data/DAList_limpa.rds` (complete matrix), `01_Trunk/05_Set_Scores/c_data/singscore.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/set_tests.rds` (read by every figure in `Figures/`), `c_data/set_contrasts.xlsx`, `b_reports/set_contrasts_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/03_Set_Contrasts/a_script/set_contrasts.R` |

```
03_Set_Contrasts/
├── a_script/set_contrasts.R
├── b_reports/set_contrasts_report.pdf
└── c_data/
    ├── set_contrasts.xlsx
    └── set_tests.rds
```

`set_tests.rds` holds two tables. `set_tests` has one row per set, contrast and method, with
`collection`, `pathway`, `n_proteins`, `direction`, `nes` (fgsea only), `p`, `fdr`, `qvalue`,
`leading_edge` (fgsea only) and `main`. `protein_results` has one row per protein and contrast,
with `logFC`, `t`, `p`, `fdr`, `pi_score` and the gene mapping.

## Workbook: `set_contrasts.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `set_summary` | a contrast, floor first | `n_tested_fgsea`, `n_fdr_fgsea`, `n_collapsed`, `n_tested_fry`, `n_fdr_fry` |
| `significant` | a set at FDR < 0.05 in a contrast, either method | `method`, `collection`, `pathway`, `direction`, `nes`, `p`, `fdr`, `qvalue`, `main`, `leading_edge` (semicolon-separated) |
| `protein_summary` | a contrast | `n_tested`, `n_up`, `n_down` (BH < 0.05), `min_fdr` |
| `fry_correlation` | the matrix fry used | `matrix`, `correlation` |

The full table, every set in every contrast, is `set_tests` in `set_tests.rds`.

## Report: `set_contrasts_report.pdf`

1. Two pages per collection. First, singscore for every set across the 45 biopsies, centred per
   set, grouped by group and time point, with fCSA mixed, mCSA, leg press 1RM and
   comp_hypertrophy above the columns (blank at T3, where no trait was measured). Second, fgsea
   NES in the nine contrasts beside a p/FDR status map. Rows are ordered within collection by
   clustering on singscore, for display only.
2. One dot plot per contrast of the fgsea sets that survived `collapsePathways()` at FDR < 0.05,
   placed by NES, coloured by collection, sized by proteins; 40 sets to a page.
3. One dot plot per contrast of every set at fry p < 0.05, placed by signed -log10 p.
