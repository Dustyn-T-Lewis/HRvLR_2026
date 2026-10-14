# 04_Continuous_no_group

Which proteins, gene sets and modules follow muscle size and strength when every person is pooled
and each trait stays a number? No group label enters any model. Traits are used as measured
(`phenotype_tables()` in `00_Setup/Shared_functions/analysis_reporting_helpers.R`): fibre area
(mixed, type I, type II), muscle area and two 1RMs as levels at T1 and T2 and as % change, plus the
composite hypertrophy score and total volume load. T3 has no phenotype.

| Walkthrough | Does | Writes |
|---|---|---|
| `01_Trait_Associations/a_script/trait_associations.qmd` | sample model, change score, two-predictor and interaction fits at three levels | `01_Proteins/`, `02_Sets/`, `03_Modules/` workbooks; `associate.rds` |
| `02_Trait_Pathways/a_script/trait_pathways.qmd` | fgsea on each trait's protein ranking | `trait_pathways.xlsx`, `pathways.rds` |
| `03_Profile_Clusters/a_script/profile_clusters.qmd` | proteins grouped by their pattern of trait links | `profile_clusters.xlsx` |

`02` and `03` read the protein change scores `01` saves, so `01` runs first.

## Inputs

`01_Processing` outputs (`DAList_corrected.rds`, `DAList_limpa.rds`, `singscore.rds`,
`gene_sets.rds`, `modules.rds`) and `00_Setup/Input/phenotype.csv`.

## Run

```sh
Rscript run_all.R 04_Continuous_no_group
```

## Outputs

Each step's `b_reports/` holds its rendered walkthrough and `c_data/` its workbooks. Every test
rests on 16 people or fewer, and HR and LR are the extremes of 25, so correlations across these 16
run stronger than they would across all 25.
