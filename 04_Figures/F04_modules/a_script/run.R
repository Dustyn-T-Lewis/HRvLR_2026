#!/usr/bin/env Rscript
# F04 WGCNA module atlas: render the two main panels (A module-trait atlas,
# B module-level NES scatters), the WGCNA-QC and phenotype supplements, then
# stitch the composite and write the workbook.
pacman::p_load(here)
source(here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))

unlink(
  setdiff(list.files(RPT_DIR, full.names = TRUE), file.path(RPT_DIR, ".gitkeep")),
  recursive = TRUE
)
dir.create(file.path(RPT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)

panel_dir <- here("04_Figures", "F04", "a_script", "panels")

# ORA first: it writes module_ora.csv, which panels A and B read for labels.
source(file.path(panel_dir, "panel_B.R"))
source(file.path(panel_dir, "panel_A.R"))
source(file.path(panel_dir, "panel_B_nes.R"))

# Supplements
source(file.path(panel_dir, "panel_C.R"))
source(file.path(panel_dir, "panel_S_phenotype_lmm.R"))
source(file.path(panel_dir, "supp_soft_threshold.R"))
source(file.path(panel_dir, "supp_dendrogram.R"))

source(here("04_Figures", "F04", "a_script", "HRvLR_F04_composite.R"))
cat("F04 complete.\n")
