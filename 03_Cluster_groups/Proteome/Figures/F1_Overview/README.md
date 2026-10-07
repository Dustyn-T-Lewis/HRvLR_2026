# F1_Overview

Figure 1 for the proteome arc, once per clustering input: who is in C1 and C2, then the
whole-proteome picture with C1 and C2 in place of HR and LR.

## The idea

The figure is the HR/LR Figure 1 (`02_HRvLR_groups/Figures/F1_Overview`) drawn five times, once
for each clustering input (T1, T2, T3, training change, acute change). Panel a changes: the design
drawing becomes a map of that input's clustering. Panels b to f are the same builders in
`R/figure_overview.R` with the groups set to `C1` and `C2`.

Here the groups come from the proteome, so some contrasts in panel b compare the clusters on the
values that formed them. Each contrast's row label carries its validity against the input
(`validity()` in `../../R/inputs.R`, explained in [`../../README.md`](../../README.md)): `(circ.)`
for circular, `(caut.)` for caution, nothing for clean.

The script reads saved outputs only. Its one computation is the phenotype PCA, a projection with
no test. Before writing anything, `figure_overview()` checks the plotted numbers against their
sources: nine contrasts whose FDR counts match `contrast_summary`, one PERMANOVA row per test,
every selected pathway in all nine contrasts, and eight barcode lanes whose FDR matches
`fry_concordance`. A mismatch stops the script.

## Panels

| Panel | Shows | Built by | Source |
|---|---|---|---|
| a | Each person clustered at that input on PC1 and PC2 of the z-scored eigengenes, coloured C1 or C2, with HR/LR as the point shape. The title names the input; the subtitle gives what was clustered, the K the gap statistic selects, the K used, and each cluster's size and bootstrap Jaccard. | `panel_clustering()` | `01_Clustering/c_data/clustering.xlsx`: `space`, `space_variance`, `choice`, `stability`, filtered to the input |
| b | Proteins per contrast, down left and up right, each row tagged with its validity. The light bar counts p < 0.05 or Π < 0.05; the dark bar, drawn over it from zero, counts the Π < 0.05 subset. The subtitle gives the number tested and at FDR < 0.05. | `panel_protein_counts()` | `04_Set_Contrasts/<input>/c_data/set_tests.rds`, `protein_results` |
| c | Proteome PCA (left) and PCA of the six size and strength traits at T1 and T2 (right), each with its PERMANOVA: time and group for the proteome, time only for the phenotype. People are relabelled C1/C2 for that input. | `panel_pca_pair()` | `03_Protein_Contrasts/<input>/c_data/protein_contrasts.xlsx`: `pca`, `pca_variance`, `permanova`; `phenotype.csv`; `groups.rds` |
| d | Pathways at fgsea FDR < 0.05 in both clusters of a window (C1 - C2: two of the three time points; Interaction: either interaction, strongest five) from Hallmark, Reactome, KEGG and GO BP, each shown in all nine contrasts. | `panel_pathway_grid()` | `set_tests.rds`, `set_tests`; `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds` |
| e, f | Each cluster's p < 0.05 sets for one window (e: training, f: acute) placed in the other cluster's ranking, drawn as limma's barcode. | `signature_barcode()` | `set_tests.rds`; `protein_contrasts.xlsx`, `fry_concordance` |

Π is Xiao et al. (2014): p^|log2 FC|, which ranks a protein on p and fold change together.

## How to read it

- Read panel b's tags first. A circular row shows the difference the clustering made; a caution
  row may show regression to the mean; only clean rows test whether the split carries over.
- Panels c to f carry no tags but follow the same rule. The proteome PERMANOVA group test uses
  every person's mean over all time points, so it shares the clustering's time points at every
  input (caution in `03_Protein_Contrasts/c_data/protein_contrasts_all_inputs.xlsx`). A pathway or
  barcode in a window inherits that window's tag in panel b.
- The phenotype played no part in the clustering, but the phenotype PCA still has no group test;
  `08_Cluster_Phenotype` tests the traits by cluster.
- C1 and C2 are PAM's cluster numbers: they carry no direction, and C1 in one input's figure is
  unrelated to C1 in another's.
- Read the Jaccard in panel a's subtitle before the rest: 0.75 or more stable, 0.6 to 0.75
  doubtful membership, below 0.6 not to be trusted.
- Panel d's rules (overlap pruning, `*` for opposite signs, the "Proteasome:" label) are those of
  the HR/LR figure and are explained in its README.

## Files

| | |
|---|---|
| Reads | per input: `03_Cluster_groups/Proteome/03_Protein_Contrasts/<input>/c_data/protein_contrasts.xlsx`, `03_Cluster_groups/Proteome/04_Set_Contrasts/<input>/c_data/set_tests.rds`, `03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/<input>/groups.rds`; shared: `03_Cluster_groups/Proteome/01_Clustering/c_data/clustering.xlsx`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `00_Input/phenotype.csv` |
| Writes | per input: `b_reports/f1_overview_<input>.pdf` (300 x 175 mm, cairo PDF), `c_data/f1_overview_<input>.xlsx` |
| Run | `Rscript 03_Cluster_groups/Proteome/Figures/F1_Overview/a_script/f1_overview.R`, or `run_all.R` after the arc's last stage |

The outputs are not tracked in git.

```
F1_Overview/
├── a_script/f1_overview.R
├── b_reports/
│   └── f1_overview_T1.pdf  _T2  _T3  _training  _acute
└── c_data/
    └── f1_overview_T1.xlsx  _T2  _T3  _training  _acute
```

## Workbook: `f1_overview_<input>.xlsx`

`read_me` lists the sheets; `columns` defines every column.

| Sheet | Holds |
|---|---|
| `clustering` | Panel a: each person's HR/LR arm, cluster as group, cluster number, PC1 and PC2. |
| `protein_counts` | Panel b: per contrast and direction, proteins at Π < 0.05, at p < 0.05 or Π, and at FDR < 0.05. |
| `pca_proteome` | Panel c, left: PC1 and PC2 per sample. |
| `pca_phenotype` | Panel c, right: PC1 and PC2 per biopsy. |
| `permanova` | Panel c: PERMANOVA for the proteome (time, group) and the phenotype (time). |
| `pathways` | Panel d: each selected set in all nine contrasts, with NES and fgsea FDR. |
| `barcode` | Panels e, f: per window and lane, set size, share moving the same way and fry FDR. |

In `protein_counts` the validity tag is part of `label`, as drawn in panel b.
