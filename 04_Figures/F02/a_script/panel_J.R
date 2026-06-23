# F02 Panel J: HR vs LR acute-response DEP overlap (Venn) ----
# Shared vs responder-specific proteins changing with the acute bout (T3 - T2), Pi < 0.05.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")
if (!exists("venn_panel")) source("04_Figures/F02/a_script/panel_I.R")

venn_panel(
  "Acute_HR", "Acute_LR",
  "Acute Response (T3 - T2)",
  file.path(RPT_DIR, "panel_j_venn_acute")
)
cat("F02 Panel J done.\n")
