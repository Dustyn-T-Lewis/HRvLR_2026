# Global ordination and PERMANOVA on all samples, strata = subject.
source(here::here("functions", "assoc_leaf.R"))
run_global_ordination_leaf(
  assoc_load(),
  here::here("04_Figures", "F04_association", "global", "ordination")
)
