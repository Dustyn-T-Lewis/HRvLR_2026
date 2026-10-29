# 04_Elastic_Net

Which proteins, modules or Hallmark sets, taken together, carry muscle size and strength, all 16
people pooled? Elastic-net regression (glmnet, α = 0.5) in six pairings of a feature value with a
trait: T1 level against the T1 trait; T1 level against the % gain; training change against the
% gain; T2 level against the T2 trait; and the acute change (the bout in trained muscle) against
the T2 trait and against the % gain. Traits: fibre area, muscle area, leg press and leg extension
1RM, and the growth score for the gains. 81 fits (pairing × trait × three feature sets).

Two readings per fit:

- stability selection (Meinshausen and Bühlmann 2010): the share of 200 subsamples of 12 people
  that keep each feature; stable at ≥ 0.6, where the expected number of chance features is about 1
  or fewer, and recurring at 0.4 to 0.6, an exploratory lead;
- held-out Q² from repeated 4-fold cross-validation, against 1,000 shuffles of the trait; a lead at
  permutation p < 0.10.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | protein values |
| `01_Processing/06_Module_Building/c_data/modules.rds` | eigengenes |
| `01_Processing/05_Set_Scores/c_data/singscore.rds` | Hallmark set scores |
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | sets for the ORA and the Hallmark list |
| `00_Setup/Input/phenotype.csv` | traits, through `trait_windows()` in `02_HRvLR/helpers.R` |

## Run

```sh
quarto render 03_Continuous/04_Elastic_Net/a_script/elastic_net.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/elastic_net.rds` | `summary` (per fit), `selection` (features kept in ≥ 40% of subsamples), `selection_ora`, `predictions` (held-out, per person); read by S6 |
| `c_data/elastic_net.xlsx` | the same tables |
| `b_reports/elastic_net.html` | the rendered walkthrough |

## Read with care

- 9 of 81 fits are leads at permutation p < 0.10 and 2 at p < 0.05, about what chance gives (8
  and 4). The strongest pairs the acute change with the muscle-area gain: modules Q² 0.24
  (p 0.008) and proteins Q² 0.14 (p 0.019), with Hallmark sets consistent (p 0.061). Its stable
  proteins (OLA1, AKR1B1, JSRP1, CSNK2B, STOML2, RPL37A, TTC37) and recurring ones (PSMA3, GSTM2,
  STBD1, EIF4E) all rise more after the bout in people who grew more (S6).
- Over-representation among selected proteins at FDR < 0.10 is a lead only: for example, the
  training change against the leg-press gain selects proteins in mitochondrial protein import and
  inner-membrane organisation (FDR 0.046).
- With 16 people, the elastic net finds features that carry a trait jointly; it is not a model for
  new people.
