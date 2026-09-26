# 02 · Differential Expression

The protein level. Tests every protein on nine contrasts, scores how well it separates the groups
on eight tasks, and ties it to phenotype per biopsy and per subject.

| Step | Runs | Writes |
|---|---|---|
| [`01_Design`](01_Design/README.md) | design formula, nine contrasts, within-subject correlation | `design.rds` |
| [`02_Contrasts`](02_Contrasts/README.md) | proteoDA limma fit, null checks, `treat()`, fry concordance | `fit.rds` |
| [`03_Classify`](03_Classify/README.md) | AUC and Wilcoxon p on eight tasks, an ROC curve per nominal protein | workbook and PDF |
| [`04_Associate`](04_Associate/README.md) | sample-level limma model per trait, change-score Spearman | `associate.rds` |

The pathway (`03_Pathway_Enrichment`) and module (`04_Network`) levels run the same three steps in
the same slots.

## Six cell means, subject random

```r
add_design(dal, "~ 0 + group + (1 | subject)")
```

One cell mean per arm and timepoint: 45 samples, 39 residual degrees of freedom. Subject is a
random effect because HR and LR are different people; a fixed subject term would absorb every
between-arm contrast. The within-subject correlation is 0.189.

| Contrast | Role | Definition |
|---|---|---|
| `Training_HR` | descriptive | HR_T2 − HR_T1 |
| `Training_LR` | descriptive | LR_T2 − LR_T1 |
| `Acute_HR` | descriptive | HR_T3 − HR_T2 |
| `Acute_LR` | descriptive | LR_T3 − LR_T2 |
| `Baseline_HRvLR` | floor | HR_T1 − LR_T1 |
| `Trained_HRvLR` | descriptive | HR_T2 − LR_T2 |
| `Acute_HRvLR` | descriptive | HR_T3 − LR_T3 |
| `Training_Interaction` | primary | (HR_T2 − HR_T1) − (LR_T2 − LR_T1) |
| `Acute_Interaction` | secondary | (HR_T3 − HR_T2) − (LR_T3 − LR_T2) |

The floor is not a negative control: the HR/LR label was cut from the training outcome, so a real
baseline difference would be predictive.

## No contrast has a BH hit

BH runs within each contrast. A protein with an empty cell in a contrast is not tested there, so
the denominators run from 1,877 to 1,892. The lowest FDR across all nine is 0.071 (Acute_HR). The
primary contrast could detect a median effect of 1.04 log2 units at 80% power.

| Contrast | prop. true null | median \|t\| | p < 0.2 |
|---|---:|---:|---:|
| Training_Interaction *(primary)* | 1.00 | 0.64 | 0.16 |
| Acute_Interaction *(secondary)* | 0.96 | 0.73 | 0.24 |
| Baseline_HRvLR *(floor)* | 1.00 | 0.61 | 0.16 |
| Training_HR | 0.90 | 0.79 | 0.25 |
| Training_LR | 1.00 | 0.64 | 0.18 |
| Acute_HR | 0.85 | 0.85 | 0.30 |
| Acute_LR | 0.98 | 0.69 | 0.24 |
| Trained_HRvLR | 1.00 | 0.63 | 0.17 |
| Acute_HRvLR | 0.94 | 0.76 | 0.24 |

A calibrated null sits at 1, 0.674 and 0.200. Five contrasts depart from it on all three:
Training_HR, Acute_HR, Acute_LR, Acute_HRvLR and Acute_Interaction. `treat()` at 1.15-fold returns
nothing. Two of the sixteen fry concordance tests reach FDR < 0.05.

## Only the acute bout clears chance as a classifier

No task has a BH hit. The HR tasks have fewer pairs than LR, so their smallest attainable p is
larger: 0.031 for HR training (6 pairs), 0.016 for HR acute (7), 0.0078 for LR (8).

| Task | Nominal / tested | Chance ratio | Smallest p |
|---|---:|---:|---:|
| Training, HR (T1 to T2) | 43 / 1,555 | 0.55 | 0.031 |
| Training, LR (T1 to T2) | 78 / 1,771 | 0.88 | 0.0078 |
| Acute bout, HR (T2 to T3) | 116 / 1,703 | 1.36 | 0.016 |
| Acute bout, LR (T2 to T3) | 107 / 1,742 | 1.23 | 0.0078 |
| HR vs LR at T1 (floor) | 37 / 1,733 | 0.43 | 0.0003 |
| HR vs LR at T2 | 61 / 1,835 | 0.66 | 0.0003 |
| HR vs LR, training change | 39 / 1,539 | 0.51 | 0.0007 |
| HR vs LR, acute change | 88 / 1,645 | 1.07 | 0.0003 |

With 6 to 8 pairs a Wilcoxon p is coarse, so fewer than 5% of null tests reach 0.05. A ratio below
1 is not worse than chance.

## Slow-fibre proteins track each person's type I share

The sample model fits every T1 and T2 biopsy (29 or 30 per trait) with its own measured trait.
Within-person, the proteins that rise most where a person's type I fibre share rose include the
slow isoforms: TPM3 ranks 2nd, TNNC1 3rd, TNNT1 5th, TNNI1 18th and MYL3 26th of 1,640. None
reaches BH < 0.05. This is the positive control: the proteome tracks a fibre-type count made
independently of the proteomics.

No model term has a BH hit. The between-person mCSA term returns 206 nominal proteins where chance
predicts 82 (ratio 2.51), the largest excess at this level; leg extension 1RM follows at 1.51
within and 1.61 between. Every other term sits between 0.56 and 1.30.

## Two change-score pairs survive BH

| Window | Outcome | Nominal / chance | BH < 0.05 |
|---|---|---:|---:|
| training | d_leg_ext_1rm | 1.96 | 0 |
| training | d_mcsa | 0.95 | 1 (RPS4X, rho 0.91, n 14) |
| baseline | d_leg_ext_1rm | 2.36 | 0 |
| baseline | d_leg_press_1rm | 1.25 | 1 (ACSL3, rho −0.87, n 15) |
| baseline | pct_leg_press_1rm | 1.63 | 0 |
| acute | comp_hypertrophy | 1.86 | 0 |
| acute | pct_mcsa | 1.82 | 0 |
| acute | pct_type1_share | 1.75 | 0 |

The 60 cells run from 0.44 to 2.36. Within-arm correlations rest on 4 to 8 subjects and carry an
exact p.
