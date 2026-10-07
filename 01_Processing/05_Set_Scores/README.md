# 05_Set_Scores

How strongly does each biopsy express each gene set?

## The idea

`04_Gene_Sets` fixed which proteins belong together. This stage turns each set into one number per
biopsy, so a set can be followed over time, compared between people, or correlated with a trait
just like a single protein.

singscore (Foroutan et al. 2018) does this with ranks. Within one sample it ranks every measured
gene from least to most abundant. A set's score is the mean rank of its members, rescaled so the
lowest and highest possible means map to the ends of the scale, and centred so that 0 means the
members sit, on average, in the middle of that sample's ranking. A positive score means the members
rank high in that sample; a negative score, low.

Three properties follow from using ranks within a sample:

- A score depends on that sample alone. Adding or dropping a biopsy changes no other biopsy's
  score.
- A score is relative. A set scores higher when its members rise against the rest of the proteome
  in that sample, not when the whole sample is brighter.
- Ranks need a value for every gene in every sample and cannot carry a covariate, so the input is
  limpa's complete matrix with the blood slope already removed (`DAList_corrected.rds`).

singscore also returns a dispersion per set and sample: the median absolute deviation of the
members' ranks. Low dispersion means the members sit together in the ranking; high means they
scatter, and the mean rank summarises them less well.

No p-value is computed here and no contrast is seen. The arcs test the scores.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Gene matrix | `filter(protein_map, selected)` | One representative protein per gene symbol from `04_Gene_Sets`; rows renamed to symbols. |
| Rank | `singscore::rankGenes()` | Ranks within each sample; ties take the lowest rank (the package default). |
| Score | `singscore::multiScore(upSetColc = sets)` | One centred score and one dispersion per qualifying set and sample. Every set is scored as an up-set: its score rises when its members rank higher. |
| Summaries | `apply()`, `rowMeans()` | Range and spread of the score matrix; mean score and dispersion per set; medians per collection. |
| Structure | `prcomp()`, `aov(~ subject)` | PCA of the samples on their set scores; the share of PC1 and PC2 variance that subject identity explains. |

## How to read the outputs

- `structure_check` shows how far the leading patterns in the scores follow the person rather than
  the time point. A large subject share is why the arcs work with within-person change beside the
  baseline level.
- `mean_sd_across_sets` is spread between sets within a sample, which reflects which genes each set
  holds. `mean_sd_across_samples` is spread of one set between samples, which is what the arcs work
  with.
- Sets that share many members get similar scores. Highly overlapping sets are not removed in
  `04_Gene_Sets`, so a result in one often appears in its neighbours.
- The figure by study group pools every set; it is a check of the score distribution, not a test.

## Files

| | |
|---|---|
| Reads | `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/03_Quantification/c_data/DAList_corrected.rds` |
| Writes | `c_data/singscore.rds` (`scores`: set by sample matrix, rows named by `set_id`), `c_data/set_scores.xlsx`, `b_reports/set_scores_report.pdf` |
| Run | `Rscript 01_Processing/05_Set_Scores/a_script/set_scores.R` |

`singscore.rds` is read by the set contrasts, set classification and trait associations in
`02_HRvLR_groups`, by `04_Continuous_no_group/01_Trait_Associations/02_Sets`, and is copied by each
`03_Cluster_groups/*/01_Clustering`.

```
05_Set_Scores/
├── a_script/set_scores.R
├── b_reports/set_scores_report.pdf
└── c_data/
    ├── set_scores.xlsx
    └── singscore.rds
```

## Workbook: `set_scores.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close it.

| Sheet | One row is | Key columns |
|---|---|---|
| `score_summary` | the whole matrix | `n_sets`, `n_samples`, `min`, `max`, `mean`, `mean_sd_across_sets`, `mean_sd_across_samples` |
| `collection_spread` | a collection | `n_sets`, `median_score`, `median_dispersion` |
| `set_spread` | a set | `collection`, `score` (mean across samples), `dispersion` (mean across samples) |
| `structure_check` | a principal component (PC1 to PC4) | `variance_explained`, `subject_share` (PC1 and PC2 only) |
| `set_scores` | a set | `set_id`, then the score in each sample |

## Report: `set_scores_report.pdf`

1. Mean score against mean dispersion, one point per set, coloured by collection.
2. The distribution of every score in every sample, pooled within each of the six arm-by-time
   groups.
