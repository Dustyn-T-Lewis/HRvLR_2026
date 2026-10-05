# 05_Associate

Ties each gene set's singscore to phenotype: in the sample model, as continuous traits with every subject pooled,
within each arm, and as differences between the arms. The code is the shared `associate_step()` in
`R/associate.R`, which the protein, set and module levels all call.

| | |
|---|---|
| Reads | `02_Scores/c_data/singscore.rds`, `01_Gene_Sets/c_data/gene_sets.rds`, `01_Preprocess/03_Quantify/c_data/DAList_corrected.rds` (sample sheet), `00_Input/phenotype.csv` |
| Writes | `c_data/05_associate.xlsx`, `b_reports/05_associate_figures.pdf` |
| Run | `Rscript 03_Pathways/05_Associate/a_script/05_associate.R` |

## Method

Traits: fibre area (mixed, type I, type II), mCSA and the leg-press and leg-extension 1RMs, as
levels at each biopsy and as % change; plus `comp_hypertrophy` and `volume_load_total_kg`. Arm is
the median split of `comp_hypertrophy`, so links with it and with fibre-area change partly carry
the group split.

Sample model. Per trait, limma fits `value ~ timepoint + between + within` on the T1 and T2 biopsies
with subject blocked: `between` is the subject's mean trait, `within` the biopsy's deviation from it.
BH within trait and term.

Change score, every subject pooled. One value per subject in five pairings: training change, T1
level and acute change against the % changes, composite and volume (`baseline` is the T1 level
against the change to come); T1 and T2 levels against the trait levels at the same biopsy
(`baseline_level`, `post_level`). Per set: limma on one value per subject (`p`, the primary
test), Spearman (`rho`, `p_rank`, from `cor.test`), biweight midcorrelation (`bicor`, WGCNA), and
the largest limma p over leave-one-subject-out refits; `robust` marks rows where all four sit below
0.05. BH and Storey q within pairing and outcome. The `chance` sheet judges each screen against 200
shuffles of whole phenotype rows across subjects: the count at p < 0.05 and the strongest hit, with
the empirical FDR (shuffle median over observed count).

Within each arm, Spearman with the exact p (`cor.test`, nine subjects or fewer) in the three change
pairings. Arm difference: `value ~ arm * outcome`, HR's slope minus LR's, with five subjects per arm;
`comp_hypertrophy` stays out because arm is its split. Two predictor: `value ~ comp_hypertrophy +
strength gain`, once per 1RM. Interaction: `value ~ growth * strength gain`, the continuous form of
the arm difference, with the leave-one-out maximum p. Flips: DGCA's Fisher z difference between the
HR and LR Spearman correlations in the training, acute, T1-level and T2-level pairings, p from 200
arm-label permutations pooled over features; `z_diff` is blank where one arm's rho is +1 or -1.

Each collection is its own screen: BH, Storey q and the shuffle counts run within collection.

## Figures

Overview pages: within- and between-person t for every trait beside status maps, then pooled rho
and limma status for every outcome, one page per pairing. Then every model hit, change-score pair,
arm-difference and two-predictor hit at p < 0.05, twelve panels to a page, ordered by p.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 per analysis and comparison. |
| `model` | Sample model: effect per trait unit, t, p, BH, Storey q per set, trait and term. |
| `change_score` | Pooled: limma effect, t, p, BH, Storey q; Spearman rho and p_rank; bicor; leave-one-out maximum p; robust. |
| `chance` | Per pairing and outcome: tested, p < 0.05, BH < 0.05, robust; 200-shuffle median and 95th percentile; empirical FDR; permutation p for the count and for the strongest hit. |
| `change_score_by_arm` | Spearman within each arm, exact p. |
| `arm_difference` | HR minus LR slope, t, p, BH, Storey q. |
| `two_predictor` | Growth and strength terms, t, p, BH, Storey q. |
| `interaction` | Growth x strength term, t, p, BH, Storey q, leave-one-out maximum p. |
| `flips` | HR and LR rho, Fisher z difference, permutation p, BH, Storey q, DGCA class. |
| `plot_index` | PDF page of each overview, and page and panel of every p < 0.05 row. |
