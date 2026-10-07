# 02_Trait_Pathways

Which gene sets lead each trait's protein ranking, with every person pooled?

| | |
|---|---|
| Reads | `04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/associate.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` |
| Writes | `c_data/pathways.rds`, `c_data/trait_pathways.xlsx`, `b_reports/trait_pathways_report.pdf` |
| Run | `Rscript 04_Continuous_no_group/02_Trait_Pathways/a_script/trait_pathways.R` |

## Key result

- 1,378 sets ranked against each of 42 trait and window pairs.
- Sets at FDR < 0.05, summed over traits: `acute` 825, `baseline` 707, `baseline_level` 963, `post_bout_level` 393, `post_level` 499, `training` 508.
- The most for one trait: `pct_leg_ext_1rm` over `baseline`, 281 sets.

## Method

- fgsea on the limma t from `01_Trait_Associations/01_Proteins` (`change_score`), one ranking per window and trait, seeded; BH within window and trait.

## Read with care

- fgsea counts run high (proteins are not exchangeable) and sets overlap, so several hits can be one signal. No self-contained test (fry) runs here.

## Sheets

`trait_pathways.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `counts` | Per window and trait: sets tested, at p < 0.05 and at FDR < 0.05 (up, down). |
| `pathways` | Per window, trait and set: fgsea NES, p and BH. |
