# 03_Set_Contrasts

Tests every set on every contrast with fgsea and fry, and marks non-redundant fgsea hits.

| | |
|---|---|
| Reads | `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `01_Trunk/05_Set_Scores/c_data/singscore.rds`, `00_Input/phenotype.csv`, `02_HRvLR_groups/02_Protein_Contrasts/c_data/fit.rds`, `02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds`, `01_Trunk/03_Quantification/c_data/DAList_limpa.rds` |
| Writes | `c_data/set_tests.rds`, `c_data/03_contrasts.xlsx`, `b_reports/03_contrasts_figures.pdf` |
| Run | `Rscript 02_HRvLR_groups/03_Set_Contrasts/a_script/03_pathways.R` |

## Method

Two methods run on every set, and they ask different questions. fgsea is competitive: does a set
rank above or below the other proteins? fry is self-contained: did the set move away from zero
across these biopsies? They can disagree when the whole proteome shifts, and fry needs more
people behind a shift than fgsea does.

**fgsea** ranks proteins by moderated t from the protein fit and asks whether a set's members sit
near one end of that list (seeded). It treats each protein as an independent draw, so a large,
tightly co-regulated set can look more significant than the number of people behind it supports.

**fry** (limma's fast rotation test) asks whether the set's proteins moved together, in one
direction, more than the model's own noise produces across these biopsies. It works on limpa's
complete matrix under the same design, subject block and within-subject correlation, and uses the
samples as the unit of evidence: 150 OXPHOS proteins rising in the same people count as one
consistent signal, not 150 confirmations. Its p therefore reflects how many subjects carry the
shift, and with seven or eight per arm that bar is high. It reports a direction: Up, Down, or Mixed
when members move both ways.

Each method's BH runs within contrast over the five collections pooled. `collapsePathways()` marks
non-redundant fgsea hits in `main` and deletes nothing. `set_tests.rds` holds `set_tests` (one row
per set, contrast and method) and `protein_results`.

## Figures

Two overview pages per collection: singscore across the 45 biopsies, centred per set, with arm, time
point, fCSA mixed, mCSA, leg press 1RM and comp_hypertrophy above the columns (the
biopsy traits are blank at T3, where none were measured); then fgsea NES in the nine contrasts
beside a status map. Then a dot plot per contrast of every fgsea collapse survivor at FDR < 0.05,
placed by NES, and one per contrast of every set at fry p < 0.05, placed by signed -log10 p, forty
sets to a page.

## Workbook

| Sheet | Holds |
|---|---|
| `set_summary` | Sets tested and called per contrast by method, and fgsea collapse survivors. |
| `significant` | Every set at FDR < 0.05, either method. The full table is `set_tests` in `set_tests.rds`. |
| `protein_summary` | Proteins tested and at BH < 0.05 per contrast. |
| `fry_correlation` | Within-subject correlation fry used, estimated on limpa's matrix. |
