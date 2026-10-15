# 03_Profile_Clusters

Groups proteins by how they relate to the phenotype. Profiles: each protein's limma t across all 42
window-trait screens, clustered (average linkage on 1 − r, dynamicTreeCut). Per trait: the p < 0.05
proteins of each training and T1-level trait, clustered on their values across people. Each group
is described by size, mean t, main module and its most over-represented gene set.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | blood-corrected protein values |
| `03_Continuous/01_Trait_Associations/c_data/associations.rds` | protein change scores |
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | sets for the over-representation test |
| `01_Processing/06_Module_Building/c_data/modules.rds` | each protein's module |
| `02_HRvLR/helpers.R` | `person_values()` |

## Run

```sh
quarto render 03_Continuous/03_Profile_Clusters/a_script/profile_clusters.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/profile_clusters.rds` | `profile_groups`, `profile_members`, `profile_means`, `trait_groups` |
| `c_data/profile_clusters.xlsx` | the same tables |
| `b_reports/profile_clusters.html` | the rendered walkthrough |

## Read with care

- 50 profile groups hold all 2,059 proteins with a mean silhouette of 0.03: the profiles have
  almost no group structure, and the groups are a way to read them, not a finding.
- 77 co-variation groups over 14 trait-window pairs, silhouettes 0.13 to 0.30.
- The over-represented sets describe groups drawn from the same t values; they do not test the
  grouping.
