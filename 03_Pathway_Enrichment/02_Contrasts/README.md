# 02_Contrasts

Tests every set on every contrast with fgsea and fry, and marks non-redundant fgsea hits.

| | |
|---|---|
| Reads | `00_Gene_Sets/c_data/gene_sets.rds`, `01_Scores/c_data/singscore.rds`, `00_Input/phenotype.csv`, `02_Differential_Expression/02_Contrasts/c_data/fit.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `c_data/set_tests.rds`, `c_data/02_contrasts.xlsx`, `b_reports/02_contrasts_figures.pdf` |
| Run | `Rscript 03_Pathway_Enrichment/02_Contrasts/a_script/02_contrasts.R` |

## Method

Both methods run on every set and both are reported; they answer different questions.

**fgsea** ranks proteins by moderated t from the protein fit and asks whether a set's members sit
near one end of that list (seeded). It treats each protein as an independent draw, so a large,
tightly co-regulated set can look highly significant when a few subjects drive it. That is why it
calls about sixty sets on `Baseline_HRvLR`, where no training difference exists.

**fry** (limma's fast rotation test) asks whether the set's proteins moved together, in one
direction, more than the model's own noise produces across these biopsies. It works on the imputed
matrix under the same design, subject block and within-subject correlation, and uses the samples as
the unit of evidence: 150 OXPHOS proteins rising in the same people count as one consistent signal,
not 150 confirmations. Its p therefore reflects how many subjects carry the shift, and with seven or
eight per arm that bar is high. It reports a direction: Up, Down, or Mixed when members move both
ways.

Read together, fgsea gives the direction and concentration of a set's change and fry says whether
enough subjects carry it. Each method's BH runs within contrast over the five collections
pooled. `collapsePathways()` marks non-redundant fgsea hits in `main` and deletes nothing.
`set_tests.rds` holds `set_tests` (one row per set, contrast and method) and `protein_results`.

ANCOVA on singscore sits beside the interactions: per set and window, `lm(later ~ earlier + arm)`
over subjects with both biopsies, T2 on T1 for training and T3 on T2 for the acute bout. BH within
window and collection.

Leave one subject out reruns the tests that feed the candidate tables once per subject dropped,
through `leave_one_out()` in `R/helpers.R`. A hit that stays at p < 0.05 in every rerun does not
rest on one person.

## Figures

Two overview pages per collection: singscore across the 45 biopsies, centred per set, with arm, time
point, fCSA mixed, mCSA, leg press 1RM and comp_hypertrophy above the columns (the
biopsy traits are blank at T3, where none were measured); then fgsea NES in the nine contrasts
beside a status map. Then a dot plot per contrast of every fgsea collapse survivor at FDR < 0.05,
placed by NES, and one per contrast of every set at fry p < 0.05, placed by signed -log10 p, forty
sets to a page. Then the ANCOVA volcano for both windows.

## Workbook

| Sheet | Holds |
|---|---|
| `set_summary` | Sets tested and called per contrast, by method, and fgsea survivors of collapsePathways. |
| `significant` | Every set at FDR < 0.05, either method. The full table is set_tests in set_tests.rds. |
| `protein_summary` | Proteins tested and at BH < 0.05 per contrast. |
| `fry_correlation` | Within-subject correlation fry used, estimated on the imputed matrix. |
| `ancova` | ANCOVA per set and window on singscore: arm effect at equal starting value, t, p, BH FDR, Storey q. |
| `leave_one_out` | Per set and HR-against-LR test: reruns with one subject left out, reruns still at p < 0.05, the largest p and the subject whose removal gave it. |
