# F02 Panel J: HR vs LR acute-response overlap, baseline-anchored (3-set Euler)
# Shared vs responder-specific proteins changing with the acute bout (T3 - T2,
# Pi < 0.05), overlaid on the baseline HR-vs-LR difference set.

pacman::p_load(here)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))
if (!exists("venn3_panel")) source(here("04_Figures", "F02", "a_script", "functions", "venn.R"))

venn3_panel(
  "Acute_HR", "Acute_LR", "Baseline_HRvLR",
  "Acute Response (T3 - T2)",
  file.path(RPT_DIR, "supp", "panel_j_venn_acute")
)
cat("F02 Panel J done.\n")
