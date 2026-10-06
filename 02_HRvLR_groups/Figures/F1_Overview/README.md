# F01_Overview

Figure 1: the design and the whole-proteome picture, before any single result.

| | |
|---|---|
| Reads | `02_HRvLR_groups/02_DEP/c_data/02_contrasts.xlsx` (`contrast_summary`, `pca`, `pca_variance`, `permanova`, `fry_concordance`), `02_HRvLR_groups/03_Pathways/c_data/set_tests.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `00_Input/phenotype.csv`, `assets/design.png` |
| Writes | `b_reports/f01_overview.pdf`, `b_reports/f01_overview.png` (300 x 175 mm, 600 dpi), `c_data/f01_overview.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/F01_Overview/a_script/f01_overview.R` |

The script reads saved outputs only. Its one computation is the phenotype PCA, a projection with
no test. Before writing anything it checks the plotted numbers against their sources: FDR counts
against `contrast_summary`, sample counts, nine columns per pathway row, and every barcode FDR
against `fry_concordance`.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Design, the measurements at each time point and the four contrast families, drawn in BioRender. | `assets/design.png` |
| b | Proteins per contrast, down left and up right. The light bar counts p < 0.05; the dark bar, drawn over it from zero, counts the Π < 0.05 subset. No protein reaches FDR < 0.05. | `set_tests.rds`, `protein_results` |
| c | Proteome PCA (left) and PCA of the six size and strength traits at T1 and T2 (right), each with a shaded 95% ellipse per time point and its PERMANOVA. The phenotype gets no arm test, because HR and LR were split on how much these traits changed. | `02_Contrasts`: `pca`, `pca_variance`, `permanova`; `phenotype.csv` |
| d | Pathways at fgsea FDR < 0.05 in both arms of a window (HR - LR: two of the three time points; Interaction: either interaction, strongest five) from Hallmark, Reactome, KEGG and GO BP, each shown in all nine contrasts. The square beside each name marks its database. | `set_tests.rds`, `set_tests` |
| e, f | Each arm's p < 0.05 sets for one window (e: training, f: acute) placed in the other arm's ranking for that window, drawn as limma's barcode (Wu & Smyth 2012) as in YvO F4F. Ticks mark members, the up set above the axis and the down set below; the curve is a tricube moving average of membership along the rank. | ranks from `set_tests.rds`; `02_Contrasts`: `fry_concordance` |

Π is Xiao et al. (2014, Bioinformatics 30:801): p^|log2 FC|. It ranks a protein on p and fold
change together and is the nominal threshold, since no protein reaches FDR < 0.05.

Panel d keeps one set from each overlapping group: sets are ranked by total |NES| across the
window's contrasts, and a set sharing half or more of its measured proteins with a stronger one is
dropped. `*` marks a set whose sign differs among the contrasts where it is significant. A set
whose significant leading edges are at least half proteasome subunits (PSM*) is labelled
"Proteasome:", because several Reactome signalling sets carry only that signal. Dot size is fgsea
FDR; fry, which rotates samples and so keeps the correlation between proteins, gives none of these
sets FDR < 0.05 (`03_Contrasts`, `set_summary`). fgsea calls 64 training and 206 acute interaction
sets, so the Interaction group is capped at five.

In e and f, the two arms are different people, so no comparison shares a group mean. Each lane
reports set size, the share of members moving the same way and fry FDR:

| Sets | Ranked in | Up | Down |
|---|---|---|---|
| LR training | HR training | 82% same way, FDR 0.016 | 76%, FDR 0.026 |
| HR training | LR training | 90%, FDR 0.14 | 74%, FDR 0.25 |
| LR acute | HR acute | 76%, FDR 0.0088 | 86%, FDR 0.0014 |
| HR acute | LR acute | 57%, FDR 0.44 | 77%, FDR 0.42 |

LR's response carries into HR's in both windows; HR's carries into LR's in the same direction
without reaching FDR 0.05. A time-linked technical effect would also move both arms together, and
the design cannot separate it from shared biology. `S03_Signatures` repeats e and f with the
Π < 0.05 sets and follows the training sets through the bout.

Colours come from `R/figure_theme.R`, sampled from the design panel's key: Training blue, Acute
orange, HR - LR slate, Interaction yellow. Strips and row boxes use the key's tint and border with
bold black text. Within Training and Acute, HR labels take a darker shade of the family colour and
LR the colour itself. Time points take the colour of the window that produced them: T1 grey, T2
blue, T3 orange.

To change panel a, edit it in BioRender, export a PNG to `assets/design.png` and rerun the script.

## Workbook

| Sheet | Holds |
|---|---|
| `protein_counts` | Panel b: proteins at Π < 0.05, at p < 0.05 or Π, and at FDR < 0.05, per contrast and direction. |
| `pca_proteome` | Panel c, left: PC1 and PC2 per sample. |
| `pca_phenotype` | Panel c, right: PC1 and PC2 per biopsy. |
| `permanova` | Panel c: PERMANOVA for the proteome (time, arm) and the phenotype (time). |
| `pathways` | Panel d: each selected set in all nine contrasts, with NES and fgsea FDR. |
| `barcode` | Panels e, f: per window and lane, set size, share moving the same way and fry FDR. |
