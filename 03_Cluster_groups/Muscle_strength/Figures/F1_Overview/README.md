# F1_Overview

Figure 1 for the strength grouping: who is in Hi and Lo, then the whole-proteome picture with Hi
and Lo in place of HR and LR.

## The idea

The figure is the HR/LR Figure 1 (`02_HRvLR_groups/Figures/F1_Overview`) drawn on this grouping's
outputs. Panel a changes: the design drawing becomes a map of the clustering, so the reader sees
who was grouped with whom before any proteome result. Panels b to f are the same builders in
`R/figure_overview.R` with the groups set to `Hi` and `Lo`.

The script reads saved outputs only. Its one computation is the phenotype PCA, a projection with
no test. Before writing anything, `figure_overview()` checks the plotted numbers against their
sources: nine contrasts whose FDR counts match `contrast_summary`, one PERMANOVA row per test,
every selected pathway in all nine contrasts, and eight barcode lanes whose FDR matches
`fry_concordance`. A mismatch stops the script.

## Panels

| Panel | Shows | Built by | Source |
|---|---|---|---|
| a | Each clustered person on PC1 and PC2 of the scaled strength changes, coloured Hi or Lo, with HR/LR as the point shape. The subtitle gives what was clustered, the K the gap statistic selects, the K used, and each group's size and bootstrap Jaccard. | `panel_clustering()` | `01_Clustering/c_data/clustering.xlsx`: `space`, `space_variance`, `choice`, `stability` |
| b | Proteins per contrast, down left and up right. The light bar counts p < 0.05 or Π < 0.05; the dark bar, drawn over it from zero, counts the Π < 0.05 subset. The subtitle gives the number tested and at FDR < 0.05. | `panel_protein_counts()` | `04_Set_Contrasts/c_data/set_tests.rds`, `protein_results` |
| c | Proteome PCA (left) and PCA of the six size and strength traits at T1 and T2 (right), each with its PERMANOVA: time and group for the proteome, time only for the phenotype. People are relabelled Hi/Lo and anyone without a group leaves the phenotype PCA. | `panel_pca_pair()` | `03_Protein_Contrasts/c_data/protein_contrasts.xlsx`: `pca`, `pca_variance`, `permanova`; `phenotype.csv`; `design.rds` (who is in each group) |
| d | Pathways at fgsea FDR < 0.05 in both groups of a window (Hi - Lo: two of the three time points; Interaction: either interaction, strongest five) from Hallmark, Reactome, KEGG and GO BP, each shown in all nine contrasts. | `panel_pathway_grid()` | `set_tests.rds`, `set_tests`; `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| e, f | Each group's p < 0.05 sets for one window (e: training, f: acute) placed in the other group's ranking, drawn as limma's barcode. | `signature_barcode()` | `set_tests.rds`; `protein_contrasts.xlsx`, `fry_concordance` |

Π is Xiao et al. (2014): p^|log2 FC|, which ranks a protein on p and fold change together.

## How to read it

- Hi and Lo were split on strength change, so the phenotype PCA has no group test; a strength
  difference between them is built in.
- The recorded reason for K = 2 is that both clusters are stable and the split is
  independent of HR/LR. Panel a shows HR/LR as shape beside Hi/Lo as colour, so how the split
  crosses HR/LR can be read by eye.
- Anyone without both strength traits has no group; they are absent from panel a and from every
  stage the figure reads.
- Read the Jaccard in panel a's subtitle before the rest: 0.75 or more stable, 0.6 to 0.75
  doubtful membership, below 0.6 not to be trusted.
- Panel d's rules (overlap pruning, `*` for opposite signs, the "Proteasome:" label) are those of
  the HR/LR figure and are explained in its README.

## Files

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/03_Protein_Contrasts/c_data/protein_contrasts.xlsx`, `03_Cluster_groups/Muscle_strength/04_Set_Contrasts/c_data/set_tests.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `00_Input/phenotype.csv`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/clustering.xlsx`, `03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds` |
| Writes | `b_reports/f1_overview.pdf` (300 x 175 mm, cairo PDF), `c_data/f1_overview.xlsx` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/Figures/F1_Overview/a_script/f1_overview.R`, or `run_all.R` after the folder's last step |

```
F1_Overview/
├── a_script/f1_overview.R
├── b_reports/f1_overview.pdf
└── c_data/f1_overview.xlsx
```

## Workbook: `f1_overview.xlsx`

`read_me` lists the sheets; `columns` defines every column.

| Sheet | Holds |
|---|---|
| `clustering` | Panel a: each person's HR/LR arm, group, cluster, PC1 and PC2. |
| `protein_counts` | Panel b: per contrast and direction, proteins at Π < 0.05, at p < 0.05 or Π, and at FDR < 0.05. |
| `pca_proteome` | Panel c, left: PC1 and PC2 per sample. |
| `pca_phenotype` | Panel c, right: PC1 and PC2 per biopsy. |
| `permanova` | Panel c: PERMANOVA for the proteome (time, group) and the phenotype (time). |
| `pathways` | Panel d: each selected set in all nine contrasts, with NES and fgsea FDR. |
| `barcode` | Panels e, f: per window and lane, set size, share moving the same way and fry FDR. |
