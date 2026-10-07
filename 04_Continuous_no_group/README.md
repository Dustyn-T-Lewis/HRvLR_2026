# 04_Continuous_no_group

Which proteins, gene sets and modules follow muscle size and strength when every person is pooled
and each trait stays a number?

The other folders split people into groups and compare the groups. Here all 16 people are pooled,
no group label enters any model, and each trait is used as measured (`phenotype_tables()` in
`R/helpers.R`): fibre area (mixed, type I, type II), muscle area and two 1RMs as levels at T1 and T2
and as % change, plus the composite hypertrophy score and total volume load as one value per
person. T3 has no phenotype.

## Windows: what meets what

| Window | Feature value per person | Trait | Asks |
|---|---|---|---|
| `training` | change, T1 to T2 | % change, composite, volume | Does the feature change more in people whose trait grew more? |
| `baseline` | T1 level | % change, composite, volume | Does the starting level predict the response to come? |
| `acute` | change, T2 to T3 | % change, composite, volume | Does the response to the final bout follow the training response? |
| `baseline_level` | T1 level | trait levels at T1 | Do people with more of the trait carry more of the feature? |
| `post_level` | T2 level | trait levels at T2 | The same after training. |
| `post_bout_level` | T3 level | trait levels at T2 | T3 follows the same training, so the T2 traits stand in; T3 also carries the bout. |

A change window needs both of its biopsies; a person missing one drops out of that window only.

## Steps

| Step | Asks | Method |
|---|---|---|
| [`01_Trait_Associations`](01_Trait_Associations/README.md) | which proteins, set scores and eigengenes track each trait | `associate_step()` in `R/associate.R` |
| [`02_Trait_Pathways`](02_Trait_Pathways/README.md) | which gene sets lead each trait's protein ranking | fgsea on the limma t from `01_Proteins` |
| [`03_Profile_Clusters`](03_Profile_Clusters/README.md) | which proteins share a pattern of trait links | clustering of each protein's t across the 42 screens |
| [`Figures/F1_Overview`](Figures/F1_Overview/README.md) | the folder on one page | saved outputs, plus a proteome PCA |

`02_Trait_Pathways` and `03_Profile_Clusters` read `01_Proteins/c_data/associate.rds`, so
`01_Proteins` runs first; `run_all.R` keeps that order.

## Read with care

- Every test rests on 16 people or fewer; read each slope with its interval.
- Pooling avoids a split, not confounding: a feature that tracks growth may track anything that
  travels with growth across these people.
