# Global mixed model across all samples: the one fit the group contrasts slice.
source(here::here("functions", "assoc_leaf.R"))
run_global_mixed_leaf(
  assoc_load(),
  here::here("04_Figures", "F04_association", "global", "mixed_model")
)
