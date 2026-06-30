# F01 Panel I: mCSA (Pre/Post x Group + Delta)
pacman::p_load(here, readr)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "F01", "a_script", "functions", "prepost.R"))

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

out <- prepost_delta_panel(
  "mCSA_Pre", "mCSA",
  expression("mCSA (" * mm^2 * ")"), expression(Delta ~ "mCSA (" * mm^2 * ")"),
  "I"
)
save_panel(out$plot, file.path(RPT, "panels", "panel_i_mcsa"), 150, 100)
write_csv(out$audit, file.path(DAT, "audit_panel_I_mcsa.csv"))
