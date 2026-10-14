# Shared analyses

These are complete analysis implementations reused by multiple substage entry scripts. They write the same results and reports as before; their inputs and output locations are passed by the calling substage. Source the required helpers and load packages in that entry script before calling these functions.

| File | Entry function | Purpose |
|------------------------|------------------------|------------------------|
| `run_design_contrasts.R` | `run_design_contrasts()` | Grouped design matrices and contrast definitions |
| `run_protein_contrasts.R` | `run_protein_contrasts()` | Protein-level comparisons |
| `run_set_contrasts.R` | `run_set_contrasts()` | Gene-set comparisons |
| `run_module_contrasts.R` | `run_module_contrasts()` | Module comparisons |
| `run_feature_discrimination.R` | `run_feature_discrimination()` | Shared feature discrimination; run_protein_discrimination(), run_set_discrimination() and run_module_discrimination() prepare each level |
| `run_continuous_associations.R` | `run_continuous_associations()` | Continuous trait associations across protein, set and module levels |
| `run_group_trait_associations.R` | `run_group_trait_associations()` | Within-group and between-group trait associations |