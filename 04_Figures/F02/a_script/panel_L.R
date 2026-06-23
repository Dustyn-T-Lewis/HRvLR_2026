# F02 Panel L (pilot): HR vs LR acute-response overlap at nominal p < 0.05 ----
# Looser-threshold companion to panel_J (Pi < 0.05); supplementary.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")
if (!exists("venn_panel")) source("04_Figures/F02/a_script/panel_I.R")

venn_panel("Acute_HR", "Acute_LR",
  "Acute Response (T3 - T2, p < 0.05)",
  metric = "p",
  file.path(RPT_DIR, "panel_l_venn_acute_p")
)
cat("F02 Panel L done.\n")
