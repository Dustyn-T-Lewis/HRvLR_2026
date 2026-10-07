# 01_Trait_Associations / 03_Modules

Which co-expression modules track each size and strength trait, with every person pooled?

## The idea

A module is a group of proteins whose abundances rise and fall together across biopsies.
`01_Processing/06_Module_Building` finds them with WGCNA, using no label or trait, and summarises each
module as an eigengene: one number per biopsy, the first principal component of the module's
proteins, signed to follow their average abundance. Proteins left unassigned (the grey module) get
no eigengene.

Each eigengene then meets each trait through the shared association step, `associate_step()` in
`R/associate.R`: the sample model, the change score in six windows, and the growth + strength and
growth x strength models. `../README.md` explains each from the ground up.

Modules turn many correlated proteins into a few tests. A module that tracks a trait says a whole
co-moving block does, which a list of single proteins cannot show directly.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Values | `readRDS()` of `modules.rds` | `eigengenes` (one row per module, one column per biopsy) and `meta` (sample sheet), checked to share one order. |
| Labels | module colour | `feature` and `label` are both the module's WGCNA colour. |
| Row order | `module_summary` | Overview rows follow module size. |
| Associations | `associate_step()` | As in `../README.md`. BH within trait and term, or window and trait. |

## How to read the outputs

- Each screen holds one test per module, so BH costs little. Below 100 tests, Storey's q is set
  equal to BH.
- An eigengene follows its module's main axis. A subset of members can follow a trait while the
  eigengene does not; `01_Proteins` shows the members one by one.
- Modules were defined on subject-centred abundance, which removes each person's own level, and
  scored on raw abundance. The between-person term therefore reads differences between people
  along an axis found within people.
- Module identity (top set, hub proteins) is in `01_Processing/07_Module_Characterisation`.

## Files

| | |
|---|---|
| Reads | `01_Processing/06_Module_Building/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/module_trait_associations.xlsx`, `b_reports/module_trait_associations_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/01_Trait_Associations/03_Modules/a_script/module_trait_associations.R` |

```
03_Modules/
├── a_script/module_trait_associations.R
├── b_reports/module_trait_associations_report.pdf
└── c_data/module_trait_associations.xlsx
```

## Workbook: `module_trait_associations.xlsx`

`read_me` lists the sheets; `columns` defines every column. `feature` and `label` are the module
colour.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | an analysis and comparison | `analysis`, `comparison`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `model` | a module, trait and term | `trait`, `term`, `n_samples`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `fdr`, `qvalue` |
| `change_score` | a module, window and trait | `window`, `trait`, `n_subjects`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `rho`, `rho_lower`, `rho_upper`, `p_rho`, `bicor`, `p_bicor`, `fdr`, `qvalue` |
| `two_predictor` | a module, window, 1RM and coefficient | `strength_test`, `coefficient`, `trait`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `fdr`, `qvalue` |
| `interaction` | a module, window and 1RM | `strength_test`, `interaction`, `interaction_lower`, `interaction_upper`, `t`, `p`, `fdr`, `qvalue` |
| `plot_index` | a report page or a drawn hit | `analysis`, `pdf_page`, `position` |

`slope` is change in eigengene per unit of trait.

## Report: `module_trait_associations_report.pdf`

1. Sample-model overview: every module, within- and between-person t per trait beside the p and
   BH status maps.
2. One change-score overview per window (six pages): Spearman rho per trait beside the limma
   status.
3. Sample-model hits, up to 24 per trait and term, twelve panels to a page.
4. Change-score hits, up to 24 per window.
5. Two-predictor hits, up to 24 per window.
