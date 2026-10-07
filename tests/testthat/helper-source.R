suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(limma)
})
source(here::here("R", "helpers.R"))
source(here::here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

# A sample table for four people over three time points, as 01_Processing writes it.
toy_meta <- function() {
  tibble(
    sample_id = paste0(rep(c("S1", "S2", "S3", "S4"), each = 3), "_", c("T1", "T2", "T3")),
    subject = rep(c("HR_S1", "HR_S2", "LR_S3", "LR_S4"), each = 3),
    arm = rep(c("HR", "HR", "LR", "LR"), each = 3),
    timepoint = rep(c("T1", "T2", "T3"), 4)
  ) |>
    mutate(group = factor(paste(arm, timepoint, sep = "_")))
}
