# 04 · Continuous, no group

Which proteins, gene sets and modules follow muscle size and strength when every person is pooled
and each trait stays a number?

## The idea

The other arcs split people into groups (HR and LR, or clusters) and compare the groups. A split
keeps only which side of a line each person falls on. This arc keeps the number: all 16 people
are pooled, no group label enters any model, and each trait is used as measured.

The traits come from `00_Input/phenotype.csv`, one row per person, read by `phenotype_tables()`
in `R/helpers.R`:

| Trait | Measured | Used as |
|---|---|---|
| `fcsa_mixed`, `fcsa_type1`, `fcsa_type2` | fibre cross-sectional area, all fibres and by type, at T1 and T2 | level at each biopsy (µm²), and % change |
| `mcsa` | muscle cross-sectional area at T1 and T2 | level (cm²), and % change |
| `leg_press_1rm`, `leg_ext_1rm` | one-repetition maximum at T1 and T2 | level (kg), and % change |
| `comp_hypertrophy` | the source's composite hypertrophy score, one value per person | change only |
| `volume_load_total_kg` | total kilograms lifted over the programme, one value per person | change only |

% change is `100 * (T2 - T1) / T1`. The traits were measured beside the T1 and T2 biopsies; T3 has
no phenotype. HR and LR are the median split of `comp_hypertrophy`, so this arc covers the ground
the group arcs cut into halves.

The arc reads three kinds of feature, all built in `01_Processing` with no label or trait:

| Level | Feature | Value per biopsy |
|---|---|---|
| Proteins | one protein | log2 abundance from limpa's complete matrix, blood slope removed |
| Sets | one gene set | singscore: where the set's proteins sit in that biopsy's ranking of all proteins |
| Modules | one WGCNA co-expression module | eigengene: one number summarising the module's proteins |

## Windows: what is paired with what

A window says which feature value each person contributes and which trait table it meets. The
first three use the change traits (six % changes, the composite and volume); the last three use
the six trait levels.

| Window | Feature value per person | Trait | Asks |
|---|---|---|---|
| `training` | change, T1 to T2 | % change, composite, volume | Does the feature change more in people whose trait grew more? |
| `baseline` | T1 level | % change, composite, volume | Does the starting level predict the response still to come? |
| `acute` | change, T2 to T3 | % change, composite, volume | Does the response to the final bout follow the training response? |
| `baseline_level` | T1 level | trait levels at T1 | At the same biopsy, do people with more of the trait carry more of the feature? |
| `post_level` | T2 level | trait levels at T2 | The same question after training. |
| `post_bout_level` | T3 level | trait levels at T2 | T3 follows the same training, so the T2 traits stand in; T3 also carries the bout's acute effect. |

A change window needs a person biopsied at both of its time points (`subject_change()`); a person
missing one drops out of that window only.

Beside the windows, a sample model uses every T1 and T2 biopsy with its own trait values and
splits each trait into a within-person part (how the person moved between their two biopsies) and
a between-person part (where the person sits against the others). `01_Trait_Associations/README.md`
explains both.

## Stages

| Stage | Asks | Method |
|---|---|---|
| [`01_Trait_Associations`](01_Trait_Associations/README.md) | which proteins, set scores and module eigengenes track each trait | `associate_step()` in `R/associate.R`: limma sample model (within and between terms); per window, limma slope on one value per person with Spearman rho and bicor beside it; growth + strength and growth x strength models |
| [`02_Trait_Pathways`](02_Trait_Pathways/README.md) | which gene sets lead each trait's protein ranking | `fgsea` on the limma t from `01_Proteins` |
| [`03_Profile_Clusters`](03_Profile_Clusters/README.md) | which proteins share a pattern of trait links | average-linkage clustering of each protein's t across every window and trait, `dynamicTreeCut`; over-representation per group |
| [`Figures/F1_Overview`](Figures/F1_Overview/README.md) | the arc on one page | saved outputs only, plus a proteome PCA |

`run_all.R` runs the three association levels, then `02_Trait_Pathways`, then
`03_Profile_Clusters`. The last two read `01_Proteins/c_data/associate.rds`, so `01_Proteins`
must run first. `F1_Overview` is not in `run_all.R`; run it by hand after the stages.
`02_HRvLR_groups/Figures/F3_People` reads the `interaction` sheet of
`protein_trait_associations.xlsx`.

## Directory tree

```
04_Continuous_no_group/
├── README.md
├── 01_Trait_Associations/
│   ├── README.md
│   ├── 01_Proteins/
│   │   ├── a_script/protein_trait_associations.R
│   │   ├── b_reports/protein_trait_associations_report.pdf
│   │   └── c_data/
│   │       ├── protein_trait_associations.xlsx
│   │       └── associate.rds
│   ├── 02_Sets/
│   │   ├── a_script/set_trait_associations.R
│   │   ├── b_reports/set_trait_associations_report.pdf
│   │   └── c_data/set_trait_associations.xlsx
│   └── 03_Modules/
│       ├── a_script/module_trait_associations.R
│       ├── b_reports/module_trait_associations_report.pdf
│       └── c_data/module_trait_associations.xlsx
├── 02_Trait_Pathways/
│   ├── a_script/trait_pathways.R
│   ├── b_reports/trait_pathways_report.pdf
│   └── c_data/
│       ├── trait_pathways.xlsx
│       └── pathways.rds
├── 03_Profile_Clusters/
│   ├── a_script/profile_clusters.R
│   ├── b_reports/profile_clusters_report.pdf
│   └── c_data/profile_clusters.xlsx
└── Figures/
    └── F1_Overview/
        ├── a_script/f1_overview.R
        ├── b_reports/f1_overview.pdf
        └── c_data/f1_overview.xlsx
```

Shared code lives in `R/`: `associate.R` (the association step), `helpers.R` (`phenotype_tables()`,
`subject_change()`, `subject_values()`, `subject_fit()`, `spearman_by_row()`, `storey_q()`,
`write_workbook()`), `figure_theme.R` and `figure_panels.R` (the figure), and
`column_dictionary.csv` (the meaning of every workbook column).

## Reading the arc

- Every test here rests on 16 people or fewer. A slope with a wide interval says little about
  size, whatever its p.
- Pooling avoids a split, not confounding. A feature that tracks growth may track anything that
  travels with growth across these people.
- Every workbook opens on `read_me` (what each sheet holds) and `columns` (what each column
  holds, from `R/column_dictionary.csv`), and ends on `input_manifest` (each file read, with its
  md5) and `package_versions`.
