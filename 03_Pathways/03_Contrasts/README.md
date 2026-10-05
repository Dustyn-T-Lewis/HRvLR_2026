# 02_Contrasts

Tests every set on every contrast with fgsea and fry, and marks non-redundant fgsea hits.

| | |
|---|---|
| Reads | `01_Gene_Sets/c_data/gene_sets.rds`, `01_Scores/c_data/singscore.rds`, `00_Input/phenotype.csv`, `02_Proteins/02_Contrasts/c_data/fit.rds`, `02_Proteins/01_Design/c_data/design.rds`, `01_Preprocess/03_Quantify/c_data/DAList_limpa.rds` |
| Writes | `c_data/set_tests.rds`, `c_data/03_contrasts.xlsx`, `b_reports/03_contrasts_figures.pdf` |
| Run | `Rscript 03_Pathways/03_Contrasts/a_script/03_contrasts.R` |

## Method

Three methods run on every set. fgsea leads, as the exploratory test; cameraPR corroborates it;
fry is the strict, people-level test.

**fgsea** ranks proteins by moderated t from the protein fit and asks whether a set's members sit
near one end of that list (seeded). It treats each protein as an independent draw, so a large,
tightly co-regulated set can look more significant than the number of people behind it supports.

**cameraPR** (limma; Wu & Smyth 2012) ranks the same moderated t but estimates how correlated each
set's members are and widens the set's error to match. A set that fgsea calls at FDR < 0.05 and
cameraPR reaches at p < 0.05 in the same direction is flagged `camera_agrees`. cameraPR adds no
hit of its own to the candidate tables. Before this step existed, a read-only check found cameraPR
kept 92% or more of fgsea's hits nominal, always in the same direction (rank correlation 0.94 to
0.97): it trims fgsea's counts without washing out its signal. It also finds sets on
`Baseline_HRvLR`, so that contrast is a comparison (the arms may differ before training), not a
clean null.

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

ANCOVA on singscore sits beside the interactions: per set and window, `lm(later ~ earlier + arm)`
over subjects with both biopsies, T2 on T1 for training and T3 on T2 for the acute bout, variances
moderated across sets (`squeezeVar`). BH within window and collection.

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
| `set_summary` | Sets tested and called per contrast by method, fgsea collapse survivors, and fgsea hits cameraPR corroborates. |
| `significant` | Every set at FDR < 0.05, any method; `camera_agrees` marks fgsea hits cameraPR corroborates. The full table is set_tests in set_tests.rds. |
| `protein_summary` | Proteins tested and at BH < 0.05 per contrast. |
| `fry_correlation` | Within-subject correlation fry used, estimated on limpa's matrix. |
| `ancova` | ANCOVA per set and window on singscore: arm effect at equal starting value, t, p, BH FDR, Storey q. |
| `leave_one_out` | Per set and HR-against-LR test: reruns with one subject left out, reruns still at p < 0.05, the largest p and the subject whose removal gave it. |
