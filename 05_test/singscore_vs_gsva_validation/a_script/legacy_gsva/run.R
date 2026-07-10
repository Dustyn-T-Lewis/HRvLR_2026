#!/usr/bin/env Rscript
# 05_test - GSVA pathway scores + mixed models (experimental, not a manuscript
# figure). Trajectory divergence (group x time) and phenotype prediction.
pacman::p_load(here)
source(here("05_test", "a_script", "setup.R"))
source(here("05_test", "a_script", "01_gsva_scores.R"))
source(here("05_test", "a_script", "02_lmm_trajectory.R"))
source(here("05_test", "a_script", "03_lmm_phenotype.R"))
source(here("05_test", "a_script", "04_multidb.R"))
source(here("05_test", "a_script", "05_parallel_figures.R"))
cat("05_test complete.\n")
