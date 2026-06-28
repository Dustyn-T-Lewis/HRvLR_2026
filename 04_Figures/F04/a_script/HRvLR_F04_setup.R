# HRvLR_F04_setup.R - Figure 4 (Training concordance: HR vs LR).
# Compares the within-group training response of high and low responders
# (Training_HR = HR T2-T1, Training_LR = LR T2-T1): where the two groups adapt
# together and where they diverge. Reads the DEP results, the imputed matrix (for
# the fry test), and the shared F03 fgsea cache (for the NES scatter - computed
# once in F03, never here). Provides dep, da, cache, pw, the contrast pair, dirs.

pacman::p_load(here, dplyr, readr, tibble)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))
source(here("04_Figures", "functions", "concordance.R"))

C_HI <- "Training_HR"
C_LO <- "Training_LR"
C_INT <- "Training_Interaction"
HI_LEVELS <- c("HR_T1", "HR_T2")
LO_LEVELS <- c("LR_T1", "LR_T2")
LABELS <- list(
  hi = "HR", lo = "LR", phase = "Training",
  x = "HR training logFC", y = "LR training logFC",
  x_short = "HR training", y_short = "LR training"
)

RPT_DIR <- here("04_Figures", "F04", "b_reports")
DAT_DIR <- here("04_Figures", "F04", "c_data")

dep <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)
da <- readRDS(here("02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"))
cache <- read_csv(here("04_Figures", "F03", "c_data", "fgsea_cache.csv"), show_col_types = FALSE)
pw <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = FALSE, exclude_variants = TRUE
)

cat(sprintf("F04 setup: %s vs %s, %d proteins, %d cached pathway rows\n", C_HI, C_LO, nrow(dep), nrow(cache)))
