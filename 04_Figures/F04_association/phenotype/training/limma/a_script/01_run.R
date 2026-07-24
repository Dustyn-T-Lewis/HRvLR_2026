# Feature ~ continuous adaptation at the training phase (moderated limma).
source(here::here("functions", "assoc_leaf.R"))
run_phenotype_leaf(
  assoc_load(),
  phase = "training",
  out_dir = here::here("04_Figures", "F04_association", "phenotype", "training", "limma")
)
