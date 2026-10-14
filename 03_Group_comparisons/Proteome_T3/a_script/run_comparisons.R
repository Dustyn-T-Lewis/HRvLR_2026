# Compares the two Proteome_T3 groups: renders each 03_Group_comparisons/code/ step.
# Rscript 03_Group_comparisons/Proteome_T3/a_script/run_comparisons.R, from the project root.
source("00_Setup/Shared_functions/run_steps.R")
run_steps("03_Group_comparisons", "Proteome_T3", c(
  "01_design_contrasts", "02_protein_contrasts", "03_set_contrasts", "04_module_contrasts",
  "05_classification", "06_trait_associations"
))
