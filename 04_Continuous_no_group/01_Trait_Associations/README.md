# 01_Trait_Associations

Which proteins, gene-set scores and module eigengenes track each size and strength trait, with
every person pooled?

## The idea

The three levels run one function, `associate_step()` in `R/associate.R`, on different matrices.
Each level's script reads its own values (proteins, set scores or eigengenes) and hands them over
with their labels; everything below is the same at every level. No group label enters.

A feature can relate to a trait in two ways, and the step asks both.

Sample model. Every T1 and T2 biopsy carries its own trait values. For each trait, the person's
mean over their two biopsies is `between`, and each biopsy's deviation from that mean is `within`.
limma then fits

```
feature ~ timepoint + between + within,  subject blocked
```

- `within` asks whether the feature rises in the biopsy where that person's trait is higher. With
  two biopsies per person, this compares each person with themselves. The `timepoint` term takes
  the shift that training gives everyone, so `within` is left with whether people whose trait
  moved more also moved more in the feature.
- `between` asks whether people with more of the trait, on average over their biopsies, carry
  more of the feature. This compares people with each other, so anything that differs between
  people and travels with the trait can drive it.

The trait enters as its level (µm², cm², kg), so the slope is change in the feature per unit of
trait. T3 has no trait values and is not in this model.

Subject blocked means the two biopsies of one person are not treated as two independent people.
`limma::duplicateCorrelation()` estimates how alike a person's biopsies are, one consensus value
across all features, and `lmFit()` uses it.

Change score. Each person contributes one feature value per window (a change or a level, see the
arc README) and one trait value from the matching table. Three statistics then relate the two
across people:

| Statistic | What it is | Why it is here |
|---|---|---|
| limma slope | `feature ~ trait` on one value per person: change in the feature per unit of trait, with its 95% interval, moderated t and p | the main test; limma borrows information on variance across all features, which steadies estimates at 16 people |
| Spearman rho | correlation of the ranks | does not assume a straight line; one extreme person moves it less |
| bicor | biweight midcorrelation (WGCNA) | a correlation that down-weights outlying points |

BH and Storey's q are computed on the limma p only. rho and bicor are checks beside it: a limma
hit whose rho is near zero rests on a few people.

Growth and strength together. Hypertrophy and strength gain are separate outcomes, and a feature
linked to one may only look linked to the other because the two vary together across people. Two
models answer this, in the three change windows (`training`, `baseline`, `acute`), with growth
as `comp_hypertrophy` and strength as the % gain in leg press or leg extension, each 1RM in its
own model:

- Two predictor: `feature ~ growth + strength`. Each slope is the feature's link to that trait
  with the other held fixed.
- Interaction: `feature ~ growth * strength`. Only the product term is kept. In this model the
  feature's slope on strength gain is `b_strength + b_interaction * growth`, so the interaction is
  how much the strength link changes per unit of growth. A positive term: the feature tracks
  strength gain more steeply in people who grew more. A negative term: less steeply. This is the
  continuous form of asking whether the strength link differs between high and low growers,
  without cutting anyone into a group. The same term reads the other way round too: how much the
  growth link changes per unit of strength gain.

Every limma estimate carries its 95% interval from `topTable(confint = TRUE)`, named after the
estimate: `slope_lower` and `slope_upper`, `interaction_lower` and `interaction_upper`.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Traits and pairings | `phenotype_tables()`, `pairing_value()` in `R/helpers.R` | Six traits as levels at T1 and T2; % changes, `comp_hypertrophy` and `volume_load_total_kg` for the change windows; the trait levels for the level windows, T2 levels for `post_bout_level`. |
| Sample model | `limma::duplicateCorrelation()`, `lmFit(block = subject)`, `eBayes(robust = TRUE)`, `topTable(confint = TRUE)` | Per trait, on the T1 and T2 biopsies where the trait was measured. Rows for the `within` and `between` terms. |
| Per-person values | `subject_change()`, `subject_values()` in `R/helpers.R` | Change T1 to T2 and T2 to T3, and the T1, T2 and T3 levels, one column per person. |
| Change score, limma | `subject_fit()` in `R/helpers.R`: `lmFit()`, `eBayes(robust = TRUE)`, `topTable(confint = TRUE)` | `feature ~ trait` per window and trait. |
| Change score, rank | `spearman_by_row()` in `R/helpers.R`: `cor.test(method = "spearman")`, `DescTools::CorCI()` | rho, its p (`p_rho`) and a 95% interval from Fisher's z transform, approximate at this n; no interval below four people. |
| Change score, robust | `WGCNA::bicorAndPvalue()` | `bicor` and `p_bicor`. |
| Two predictor | `subject_fit(~ growth + strength)` | Both slopes, per change window and 1RM. |
| Interaction | `subject_fit(~ growth * strength)` | The `growth:strength` term, per change window and 1RM. |
| Multiple testing | `p.adjust(method = "BH")`, `storey_q()` | `fdr` and `qvalue` within each screen (below). Storey's q equals BH when a screen holds fewer than 100 tests. |
| Report and workbook | `ComplexHeatmap`, `ggplot2`, `write_workbook()` | See each level's README. |

The screens BH runs within:

| Analysis | One screen is |
|---|---|
| Sample model | one trait and term |
| Change score | one window and trait |
| Two predictor | one window, 1RM and coefficient |
| Interaction | one window and 1RM |

At the set level each collection is also its own screen.

## How to read the outputs

- Read the slope with its interval. At 16 people the intervals are wide, and a nominal p < 0.05
  says little about size.
- `within` is the stronger evidence for a link to change, because each person is their own
  control. `between` can reflect any stable difference between people.
- Each window and trait is its own screen. A feature at FDR < 0.05 in one window and not another
  does not show the two windows differ.
- `baseline` and `baseline_level` share the same T1 feature values; they differ only in the
  trait they meet. `training`, `baseline` and `acute` all meet the same change traits.
- The `acute` window uses training-change traits: it asks whether the bout response follows
  how the person responded to training, not whether the bout changed any trait.
- `post_bout_level` mixes trained state with the bout's acute effect, and its traits were
  measured at T2, not T3.
- The two 1RMs are never combined; read leg press and leg extension as separate tests.
- The interaction product rests most on the people at the extremes of both growth and strength
  gain. With 16 people, one person can carry it.
- Counts of features at p < 0.05 mean little alone: by chance, about 5% of features land there.

## Levels

| Level | Values | Rows ordered by (display only) | Differs from the others |
|---|---|---|---|
| [`01_Proteins`](01_Proteins/README.md) | limpa's complete blood-corrected matrix | average-linkage clustering of limpa's matrix | carries `observed_share`; saves `associate.rds` for `02_Trait_Pathways` and `03_Profile_Clusters` |
| [`02_Sets`](02_Sets/README.md) | singscore per qualifying gene set | the same clustering on singscore, within collection | BH and q within collection; one block of overview pages per collection |
| [`03_Modules`](03_Modules/README.md) | WGCNA module eigengenes | module size | none |

```
01_Trait_Associations/
├── README.md
├── 01_Proteins/
│   ├── a_script/protein_trait_associations.R
│   ├── b_reports/protein_trait_associations_report.pdf
│   └── c_data/
│       ├── protein_trait_associations.xlsx
│       └── associate.rds
├── 02_Sets/
│   ├── a_script/set_trait_associations.R
│   ├── b_reports/set_trait_associations_report.pdf
│   └── c_data/set_trait_associations.xlsx
└── 03_Modules/
    ├── a_script/module_trait_associations.R
    ├── b_reports/module_trait_associations_report.pdf
    └── c_data/module_trait_associations.xlsx
```

## Workbook sheets every level shares

`read_me` lists the sheets; `columns` defines every column. Sets add `collection`; proteins add
`observed_share`.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | an analysis and comparison | `analysis` (for example `model: within`, `change score: training`, `two predictor: acute`, `interaction: baseline`), `comparison`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `model` | a feature, trait and term | `trait`, `term`, `n_samples`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `fdr`, `qvalue` |
| `change_score` | a feature, window and trait | `window`, `trait`, `n_subjects`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `rho`, `rho_lower`, `rho_upper`, `p_rho`, `bicor`, `p_bicor`, `fdr`, `qvalue` |
| `two_predictor` | a feature, window, 1RM and coefficient | `window`, `strength_test`, `coefficient` (`growth` or `strength`), `trait`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `fdr`, `qvalue` |
| `interaction` | a feature, window and 1RM | `window`, `strength_test`, `interaction`, `interaction_lower`, `interaction_upper`, `t`, `p`, `fdr`, `qvalue` |
| `plot_index` | a report page or a drawn hit | `analysis`, `pdf_page`, `position`, with the window, trait and feature it draws |

Every row carries `feature` (the identifier) and `label` (a readable name).

## Report pages every level shares

1. Overviews, all features one row each. A sample-model page: within- and between-person t for
   every trait, beside their p and BH status maps. Then one page per window: Spearman rho with
   every trait, beside the limma p and BH status. Colours on the t maps stop at the 99th
   percentile of |t|.
2. Sample-model hits: every p < 0.05 row, up to 24 per trait and term, twelve panels to a page,
   smallest p first. Each panel plots the feature against the trait at the biopsy; T1 open, T2
   filled, a line joining each person's two biopsies. The within term reads those lines' slopes,
   the between term the spread of people.
3. Change-score hits: every p < 0.05 feature-trait pair, up to 24 per window, one point per
   person with a least-squares line as a guide.
4. Two-predictor hits: up to 24 per window. The x axis shows the term's own trait; the drawing
   does not hold the other predictor fixed, though the slope in the header does.

The report draws no interaction rows; read them in the workbook. The workbook keeps every row
that the 24-hit cap leaves off the pages, and `plot_index` gives the page and panel of each drawn
one.
