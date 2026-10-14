# 05_Figures

Figures drawn from saved results only. Any computation a figure adds is a description with no new
test (F1's phenotype PCA, F3's Spearman across proteins, S02's Spearman of NES across sets). Every
figure is a plain R script: it writes `b_reports/<name>.pdf` and `c_data/<name>.xlsx`, the numbers
each panel draws, and stops if a plotted number differs from its source.

| Folder | Status | Asks |
|---|---|---|
| `<grouping>/F1_Overview` | active | Who is compared, and how does the whole proteome respond before any single result? One per grouping, from `code/f1_overview.R` |
| [`HRvLR/F2_Programmes`](HRvLR/F2_Programmes/README.md) | parked draft | Which co-expression programmes move with training and the bout, and do HR and LR differ in them? |
| [`HRvLR/F3_People`](HRvLR/F3_People/README.md) | parked draft | Can one feature tell HR from LR, and does a protein's link to strength gain differ between them? |
| [`HRvLR/S01_Volcano`](HRvLR/S01_Volcano/README.md) | active | Which proteins and fgsea sets move in each within-group and interaction contrast? |
| [`HRvLR/S02_NES_Scatter`](HRvLR/S02_NES_Scatter/README.md) | active | Do HR and LR shift the same gene sets, set by set? |
| [`HRvLR/S03_Signatures`](HRvLR/S03_Signatures/README.md) | active | Supplement to F1: the signature ranks with the stricter Π sets, and the training proteins through the bout |
| [`Processing/S01_QC`](Processing/S01_QC/README.md) | active | Did the samples hold up, and does the analysis behave? |
| [`Continuous/F1_Overview`](Continuous/F1_Overview/README.md) | active | The continuous analysis on a page |

Parked means the script runs and its outputs are rebuilt, but the figure is not final. Colours,
the contrast palette and panel builders come from `00_Setup/Shared_functions/figure_style.R`,
`figure_panels.R` and `group_overview_figure.R`.

## F1, one per grouping

`code/f1_overview.R` defines `f1_overview(grouping)`; each `<grouping>/F1_Overview/a_script/` runs
it for its grouping. HRvLR's panel a is the BioRender design drawing; see
[`HRvLR/F1_Overview`](HRvLR/F1_Overview/README.md) for its panels in full. For the seven clustered
groupings panel a maps the clustering instead, and the groups are HR and LR (muscle) or C1 and C2
(proteome).

| Panel | Shows | Built by | Source, in `03_Group_comparisons/<grouping>/c_data/` unless named |
|---|---|---|---|
| a | Each clustered person on PC1 and PC2 of what they were clustered on, coloured by group, with HR/LR as the point shape. The subtitle gives what was clustered, the K the gap statistic selects, the K used, and each group's size and bootstrap Jaccard. | `panel_clustering()` | `02_Groups/<grouping>/c_data/clustering.xlsx`: `space`, `space_variance`, `choice`, `stability` |
| b | Proteins per contrast, down left and up right. The light bar counts p < 0.05 or Π < 0.05; the dark bar, drawn over it from zero, counts the Π < 0.05 subset. The subtitle gives the number tested and at FDR < 0.05. | `panel_protein_counts()` | `03_set_tests.rds`, `protein_results` |
| c | Proteome PCA (left) and PCA of the six size and strength traits at T1 and T2 (right), each with its PERMANOVA: time and group for the proteome, time only for the phenotype. | `panel_pca_pair()` | `02_protein_contrasts.xlsx`: `pca`, `pca_variance`, `permanova`; `phenotype.csv`; `01_design.rds` |
| d | Pathways at fgsea FDR < 0.05 in both groups of a window (group difference: two of the three time points; interaction: either interaction, strongest five) from Hallmark, Reactome, KEGG and GO BP, each shown in all nine contrasts. | `panel_pathway_grid()` | `03_set_tests.rds`, `set_tests`; `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| e, f | Each group's p < 0.05 proteins for one window (e: training, f: acute) placed in the other group's ranking, drawn as limma's barcode. | `signature_barcode()` | `03_set_tests.rds`; `02_protein_contrasts.xlsx`, `fry_concordance` |

`c_data/f1_overview.xlsx` holds one sheet per panel: `clustering` (a), `protein_counts` (b),
`pca_proteome`, `pca_phenotype` and `permanova` (c), `pathways` (d), `barcode` (e, f).

### Read with care

- Read the Jaccard in panel a's subtitle first: 0.75 or more stable, 0.6 to 0.75 doubtful, below
  0.6 not to be trusted. Every clustering found K = 1; K = 2 was imposed.
- Muscle_size and Muscle_strength were split on size or strength change, so the phenotype PCA has no
  group test and a difference in that trait is built in. Muscle_size reproduces the HR/LR split
  exactly, so its panels b to f repeat the HR/LR figure; Muscle_strength is stable
  and unrelated to HR/LR.
- Proteome groupings tag each contrast in panel b with its validity against the clustering:
  `(circ.)` circular, built from the same time points; `(caut.)` caution, shares a time point;
  untagged rows are clean. Only clean rows test whether the split carries over. Panels c to f carry
  no tags but follow the same rule; the proteome PERMANOVA group test uses every time point, so it is
  caution at every input.
- C1 and C2 are PAM's numbers: no direction, and unrelated between proteome groupings.
  `03_Group_comparisons/Proteome_clusters` tests the traits by cluster.

## Run

From the project root, after `03_Group_comparisons` and `04_Continuous_no_group`:

```sh
Rscript run_all.R 05_Figures                                  # every figure
Rscript 05_Figures/Proteome_T1/F1_Overview/a_script/f1_overview.R    # one figure
```
