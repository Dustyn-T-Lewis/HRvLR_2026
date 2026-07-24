# Feature ~ continuous adaptation at the baseline phase (moderated limma).
source(here::here("functions", "assoc_leaf.R"))
run_phenotype_leaf(
  assoc_load(),
  phase = "baseline",
  out_dir = here::here("04_Figures", "F04_association", "phenotype", "baseline", "limma")
)
