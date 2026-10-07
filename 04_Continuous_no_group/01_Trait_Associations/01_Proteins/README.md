# 01_Trait_Associations / 01_Proteins

Which proteins track each size and strength trait, with every person pooled?

## The idea

Each protein's log2 abundance meets each trait through the shared association step,
`associate_step()` in `R/associate.R`: the sample model (within- and between-person terms), the
change score in six windows, and the growth + strength and growth x strength models.
`../README.md` explains each from the ground up; this file covers what is particular to proteins.

The values are limpa's complete matrix with the blood slope removed (`DAList_corrected.rds`).
limpa gives every protein a value in every biopsy, filling a missing one from a detection curve
fitted across the run. Rank tests and correlations take no covariate, so this level reads the
matrix with blood already taken out rather than fitting blood as a term. The limma fits here use
those values as they are, with no precision weights, so a filled value counts as much as a
measured one. `observed_share` says how much of each protein was measured.

Gene symbols label the rows; the script stops if any protein lacks one or two share one.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Values | `readRDS()` of `DAList_corrected.rds` | `data` (log2 abundance, every biopsy), `metadata` (subject and time point), `annotation` (gene symbol), `observed` (measured or filled). |
| Measured share | `rowMeans(dal$observed)` | `observed_share` per protein. |
| Row order | `hclust(as.dist(1 - cor(t(...))), "average")` on limpa's matrix (`DAList_limpa.rds`) | Orders the overview heatmaps only; no test reads it. |
| Associations | `associate_step()` | As in `../README.md`. BH within trait and term, or window and trait. |
| Save for later stages | `saveRDS()` (`save_change_score = TRUE`) | `associate.rds` holds the `change_score` table, read by `02_Trait_Pathways` (ranks by its t) and `03_Profile_Clusters` (clusters its t profiles). |

## How to read the outputs

- Check `observed_share` before trusting a hit. A protein measured in few biopsies has most of its
  values filled, and a slope built on filled values reflects the detection curve as much as the
  person.
- A protein's link to a trait across 16 people is a correlation. It does not say the protein
  drives growth or strength.
- Proteins move together. Twenty hits from one complex are one signal seen twenty times; the
  module level (`03_Modules`) and set level (`02_Sets`) group them.
- `02_HRvLR_groups/Figures/F3_People` reads the `interaction` sheet, so a change here moves that
  figure.

## Files

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/DAList_corrected.rds` (values, annotation, sample sheet, observed mask), `01_Processing/03_Quantification/c_data/DAList_limpa.rds` (overview row order), `00_Input/phenotype.csv` |
| Writes | `c_data/protein_trait_associations.xlsx`, `c_data/associate.rds`, `b_reports/protein_trait_associations_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/01_Trait_Associations/01_Proteins/a_script/protein_trait_associations.R` |

```
01_Proteins/
├── a_script/protein_trait_associations.R
├── b_reports/protein_trait_associations_report.pdf
└── c_data/
    ├── protein_trait_associations.xlsx
    └── associate.rds
```

## Workbook: `protein_trait_associations.xlsx`

`read_me` lists the sheets; `columns` defines every column. `feature` is the UniProt accession and
`label` the gene symbol; every result sheet carries `observed_share` after `label`.

| Sheet | One row is | Key columns |
|---|---|---|
| `counts` | an analysis and comparison | `analysis`, `comparison`, `n_tested`, `n_nominal`, `n_fdr_05` |
| `model` | a protein, trait and term | `trait`, `term`, `n_samples`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `fdr`, `qvalue` |
| `change_score` | a protein, window and trait | `window`, `trait`, `n_subjects`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `rho`, `rho_lower`, `rho_upper`, `p_rho`, `bicor`, `p_bicor`, `fdr`, `qvalue` |
| `two_predictor` | a protein, window, 1RM and coefficient | `strength_test`, `coefficient`, `trait`, `slope`, `slope_lower`, `slope_upper`, `t`, `p`, `fdr`, `qvalue` |
| `interaction` | a protein, window and 1RM | `strength_test`, `interaction`, `interaction_lower`, `interaction_upper`, `t`, `p`, `fdr`, `qvalue` |
| `plot_index` | a report page or a drawn hit | `analysis`, `pdf_page`, `position` |

`associate.rds` is a list with one element, `change_score`, the same table as the sheet.

## Report: `protein_trait_associations_report.pdf`

1. Sample-model overview: every protein, within- and between-person t per trait beside the p and
   BH status maps.
2. One change-score overview per window (six pages): Spearman rho per trait beside the limma
   status.
3. Sample-model hits, up to 24 per trait and term, twelve panels to a page.
4. Change-score hits, up to 24 per window.
5. Two-predictor hits, up to 24 per window.

Panel headers give the slope or rho, p, BH q and the number of biopsies or people.
