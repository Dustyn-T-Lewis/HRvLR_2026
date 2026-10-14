# F1_Overview

Figure 1 for the strength grouping: who is in Hi and Lo, then the whole-proteome picture with Hi and Lo in place of HR and LR. Panels b to f are the HR/LR Figure 1 builders in `00_Setup/Shared_functions/group_overview_figure.R`; panel a maps the clustering in place of the design drawing.

| Panel | Shows | Built by | Source |
|---|---|---|---|
| a | Each clustered person on PC1 and PC2 of the scaled strength changes, coloured Hi or Lo, with HR/LR as the point shape. The subtitle gives what was clustered, the K the gap statistic selects, the K used, and each group's size and bootstrap Jaccard. | `panel_clustering()` | `01_Clustering/c_data/clustering.xlsx`: `space`, `space_variance`, `choice`, `stability` |
| b | Proteins per contrast, down left and up right. The light bar counts p < 0.05 or Π < 0.05; the dark bar, drawn over it from zero, counts the Π < 0.05 subset. The subtitle gives the number tested and at FDR < 0.05. | `panel_protein_counts()` | `04_Set_Contrasts/c_data/set_tests.rds`, `protein_results` |
| c | Proteome PCA (left) and PCA of the six size and strength traits at T1 and T2 (right), each with its PERMANOVA: time and group for the proteome, time only for the phenotype. People are relabelled Hi/Lo and anyone without a group leaves the phenotype PCA. | `panel_pca_pair()` | `03_Protein_Contrasts/c_data/protein_contrasts.xlsx`: `pca`, `pca_variance`, `permanova`; `phenotype.csv`; `design.rds` (who is in each group) |
| d | Pathways at fgsea FDR < 0.05 in both groups of a window (Hi - Lo: two of the three time points; Interaction: either interaction, strongest five) from Hallmark, Reactome, KEGG and GO BP, each shown in all nine contrasts. | `panel_pathway_grid()` | `set_tests.rds`, `set_tests`; `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| e, f | Each group's p < 0.05 proteins for one window (e: training, f: acute) placed in the other group's ranking, drawn as limma's barcode. | `signature_barcode()` | `set_tests.rds`; `protein_contrasts.xlsx`, `fry_concordance` |

| | |
|---|---|
| Reads | `03_Cluster_groups/Muscle_strength/03_Protein_Contrasts/c_data/protein_contrasts.xlsx`, `03_Cluster_groups/Muscle_strength/04_Set_Contrasts/c_data/set_tests.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `00_Setup/Input/phenotype.csv`, `03_Cluster_groups/Muscle_strength/01_Clustering/c_data/clustering.xlsx`, `03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds` |
| Writes | `b_reports/f1_overview.pdf` (300 x 175 mm, cairo PDF), `c_data/f1_overview.xlsx` |
| Run | `Rscript 03_Cluster_groups/Muscle_strength/Figures/F1_Overview/a_script/f1_overview.R`, or `run_all.R` after the folder's last step |

## Read with care

- Hi and Lo were split on strength change, so the phenotype PCA has no group test and a strength difference between them is built in.
- K = 2 was imposed because both clusters are stable and the split is unrelated to HR/LR. Panel a shows HR/LR as shape beside Hi/Lo as colour.
- Read the Jaccard in panel a's subtitle first: 0.75 or more stable, 0.6 to 0.75 doubtful, below 0.6 not to be trusted.
- Panel d's selection rules are those of the HR/LR figure (`02_HRvLR_groups/Figures/F1_Overview`). The script stops if a plotted count differs from its source.

## Sheets

| Sheet | Holds |
|---|---|
| `clustering` | Panel a: each person's HR/LR arm, group, cluster, PC1 and PC2. |
| `protein_counts` | Panel b: per contrast and direction, proteins at Π < 0.05, at p < 0.05 or Π, and at FDR < 0.05. |
| `pca_proteome` | Panel c, left: PC1 and PC2 per sample. |
| `pca_phenotype` | Panel c, right: PC1 and PC2 per biopsy. |
| `permanova` | Panel c: PERMANOVA for the proteome (time, group) and the phenotype (time). |
| `pathways` | Panel d: each selected set in all nine contrasts, with NES and fgsea FDR. |
| `barcode` | Panels e, f: per window and lane, set size, share moving the same way and fry FDR. |
