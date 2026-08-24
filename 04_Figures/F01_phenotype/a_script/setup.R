# F01 setup: style, the phenotype table, and the output paths. Panels and
# 01_run_phenotype.R source this first; idempotent, writes nothing.
#
# Provides: pheno, PHENOTYPES, TRAIT_LABELS, CHANGE_TRAITS, F01_RPT, F01_DAT
pacman::p_load(here, dplyr, tibble, readr)

source(here("functions", "shared_style.R"))
source(here("functions", "association.R"))

F01_RPT <- here("04_Figures", "F01_phenotype", "b_reports")
F01_DAT <- here("04_Figures", "F01_phenotype", "c_data")

pheno <- phenotype_table()
PHENOTYPES <- setdiff(names(pheno), c("subject", "group_arm"))

TRAIT_LABELS <- c(
  comp_hypertrophy = "Composite hypertrophy",
  d_fcsa_I = "fCSA type I", d_fcsa_II = "fCSA type II",
  d_fcsa_mixed = "fCSA mixed", d_nfibre_mixed = "Fibre count, mixed",
  d_nfibre_I = "Fibre count, type I", d_mcsa = "Whole-muscle CSA",
  d_1rm_legpress = "1RM leg press", d_1rm_ext = "1RM leg extension",
  volume_load = "Training volume load"
)

# volume_load is a total, not a change, so "did it move" is not a question that
# applies to it. It still belongs in the correlation structure.
CHANGE_TRAITS <- setdiff(PHENOTYPES, c("comp_hypertrophy", "volume_load"))

if (!exists("F01_PANELS")) F01_PANELS <- list()
if (!exists("F01_AUDIT")) F01_AUDIT <- list()
