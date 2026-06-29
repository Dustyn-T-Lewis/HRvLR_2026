# HRvLR_F05_setup.R - Figure 5 (Acute concordance: HR vs LR).
# Compares the acute (last-bout) response of high and low responders
# (Acute_HR = HR T3-T2, Acute_LR = LR T3-T2). Same builders as F04, Acute pair.
# Reads the DEP results, the imputed matrix (fry), and the shared F03 fgsea cache.

pacman::p_load(here, dplyr, readr, tibble)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))
source(here("04_Figures", "functions", "concordance.R"))

C_HI <- "Acute_HR"
C_LO <- "Acute_LR"
C_INT <- "Acute_Interaction"
HI_LEVELS <- c("HR_T2", "HR_T3")
LO_LEVELS <- c("LR_T2", "LR_T3")
LABELS <- list(
  hi = "HR", lo = "LR", phase = "Acute",
  x = "HR acute logFC", y = "LR acute logFC",
  x_short = "HR acute", y_short = "LR acute"
)

RPT_DIR <- here("04_Figures", "F05", "b_reports")
DAT_DIR <- here("04_Figures", "F05", "c_data")

dep <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)
da <- readRDS(here("02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"))
cache <- as_tibble(openxlsx::read.xlsx(
  here("04_Figures", "F03", "c_data", "F03_source_data.xlsx"),
  sheet = "fgsea_all"
))
pw <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = FALSE, exclude_variants = TRUE
)

cat(sprintf("F05 setup: %s vs %s, %d proteins, %d cached pathway rows\n", C_HI, C_LO, nrow(dep), nrow(cache)))
