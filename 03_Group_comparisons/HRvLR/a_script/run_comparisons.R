# Compares HR with LR: renders each 03_Group_comparisons/code/ step for this grouping.
# Rscript 03_Group_comparisons/HRvLR/a_script/run_comparisons.R, from the project root.
source("00_Setup/Shared_functions/run_steps.R")
run_steps("03_Group_comparisons", "HRvLR", c(
  "01_design_contrasts", "02_protein_contrasts", "03_set_contrasts", "04_module_contrasts",
  "05_classification", "06_trait_associations", "07_module_preservation"
))
