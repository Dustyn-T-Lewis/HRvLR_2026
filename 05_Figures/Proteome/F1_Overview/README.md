# F1_Overview

Figure 1 for the proteome grouping, once per clustering input (T1, T2, T3, training, acute): who is in C1 and C2, then the HR/LR Figure 1 builders with C1 and C2 in place of HR and LR. Each contrast in panel b carries its validity against the input: `(circ.)` circular, `(caut.)` caution, nothing for clean ([`../../README.md`](../../README.md)).

| Panel | Shows | Built by | Source |
|---|---|---|---|
| a | Each person clustered at that input on PC1 and PC2 of the z-scored eigengenes, coloured C1 or C2, with HR/LR as the point shape. The title names the input; the subtitle gives what was clustered, the K the gap statistic selects, the K used, and each cluster's size and bootstrap Jaccard. | `panel_clustering()` | `01_Clustering/c_data/clustering.xlsx`: `space`, `space_variance`, `choice`, `stability`, filtered to the input |
| b | Proteins per contrast, down left and up right, each row tagged with its validity. The light bar counts p < 0.05 or Π < 0.05; the dark bar, drawn over it from zero, counts the Π < 0.05 subset. The subtitle gives the number tested and at FDR < 0.05. | `panel_protein_counts()` | `04_Set_Contrasts/<input>/c_data/set_tests.rds`, `protein_results` |
| c | Proteome PCA (left) and PCA of the six size and strength traits at T1 and T2 (right), each with its PERMANOVA: time and group for the proteome, time only for the phenotype. People are relabelled C1/C2 for that input. | `panel_pca_pair()` | `03_Protein_Contrasts/<input>/c_data/protein_contrasts.xlsx`: `pca`, `pca_variance`, `permanova`; `phenotype.csv`; `design.rds` (who is in each group) |
| d | Pathways at fgsea FDR < 0.05 in both clusters of a window (C1 - C2: two of the three time points; Interaction: either interaction, strongest five) from Hallmark, Reactome, KEGG and GO BP, each shown in all nine contrasts. | `panel_pathway_grid()` | `set_tests.rds`, `set_tests`; `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| e, f | Each cluster's p < 0.05 proteins for one window (e: training, f: acute) placed in the other cluster's ranking, drawn as limma's barcode. | `signature_barcode()` | `set_tests.rds`; `protein_contrasts.xlsx`, `fry_concordance` |

| | |
|---|---|
| Reads | per input: `03_Group_comparisons/Proteome/02_Protein_Contrasts/<input>/c_data/protein_contrasts.xlsx`, `03_Group_comparisons/Proteome/03_Set_Contrasts/<input>/c_data/set_tests.rds`, `03_Group_comparisons/Proteome/01_Design_Contrasts/<input>/c_data/design.rds`; shared: `03_Group_comparisons/Proteome/01_Clustering/c_data/clustering.xlsx`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `00_Setup/Input/phenotype.csv` |
| Writes | per input: `b_reports/f1_overview_<input>.pdf` (300 x 175 mm, cairo PDF), `c_data/f1_overview_<input>.xlsx` |
| Run | `Rscript 05_Figures/Proteome/F1_Overview/a_script/f1_overview.R`, or `run_all.R` after the folder's last step |

The outputs are not tracked in git.

## Read with care

- Read panel b's tags first: a circular row shows the difference the clustering made, a caution row may show regression to the mean, only clean rows test whether the split carries over.
- Panels c to f carry no tags but follow the same rule; the proteome PERMANOVA group test uses every time point, so it is caution at every input.
- C1 and C2 are PAM's numbers: no direction, and unrelated between inputs. Read the Jaccard in panel a's subtitle first.
- `08_Cluster_Phenotype` tests the traits by cluster; the phenotype PCA has no group test.

## Sheets

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
