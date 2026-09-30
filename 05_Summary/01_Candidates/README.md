# 01_Candidates

Lists every protein, set and module with at least one p < 0.05 among the tests that bear on HR
against LR, one table for training and one for the acute bout.

| | |
|---|---|
| Reads | the contrasts, classify and associate workbooks of stages 02 to 04, `03_Pathway_Enrichment/02_Contrasts/c_data/set_tests.rds`, `03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds`, `04_Network/01_Build_Modules/c_data/01_build_modules.xlsx` |
| Writes | `c_data/01_candidates.xlsx` |
| Run | `Rscript 05_Summary/01_Candidates/a_script/01_candidates.R` |

## Method

Each feature collects the tests below; every p < 0.05 counts once and nothing is weighted. Each hit
carries its BH FDR and Storey q from its own screen, as exploratory corrections; `n_bh_20` and
`n_storey_20` count hits below 0.20 and `min_bh` and `min_storey` give the smallest. Rows are
sorted by hit count, then by the smallest p.

| Family | Training list | Acute list |
|---|---|---|
| Contrasts | `Training_Interaction`, `Trained_HRvLR`, `Acute_Interaction`, ANCOVA training | `Acute_Interaction`, `Acute_HRvLR`, ANCOVA acute |
| Classify | Wilcoxon on training change and at T2 | Wilcoxon on acute change |
| Associate | training window: pooled rho, rho within HR and within LR, arm slope difference, two-predictor growth and strength terms, over the eight outcomes | the same in the acute window |

Contrasts mean the limma logFC for proteins, fgsea and fry for sets, and eigengene limma and fry for
modules. `Acute_Interaction` feeds both lists. A feature on both lists is flagged `on_both`.

Hit counts overlap by construction: the interaction, the training-change classification and the
`comp_hypertrophy` association all read each subject's training change, and arm is the median split
of `comp_hypertrophy`. The `hit_distribution` sheet gives the reference: the same rules at T1
(`Baseline_HRvLR`, the T1 classification and the T1-level associations), before any training.
Baseline logFC and p, and for proteins the interaction p before blood adjustment, sit beside each row
and filter nothing. `n_stable` counts the contrast hits that stay at p < 0.05 in every
leave-one-subject-out rerun (each level's `02_Contrasts`, sheet `leave_one_out`).

## Figures

None.

## Workbook

| Sheet | Holds |
|---|---|
| `proteins_training`, `proteins_acute` | Proteins with at least one hit, by hit count then smallest p; `hits` lists each test with its effect and p. |
| `sets_training`, `sets_acute` | The same for gene sets. |
| `modules_training`, `modules_acute` | The same for modules. |
| `hit_distribution` | Features by hit count per level and list, including the baseline reference, with tests per feature. |
| `evidence` | Every p < 0.05 behind the tables, one row per test. |
