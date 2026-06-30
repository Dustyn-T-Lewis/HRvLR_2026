# F01 Panel H: 1RM Extension (Pre/Post x Group + Delta)
pacman::p_load(here, readr)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "F01", "a_script", "functions", "prepost.R"))

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

out <- prepost_delta_panel(
  "X1RM._Ext_Pre", "1RM Extension", "1RM (kg)", expression(Delta ~ "1RM (kg)"), "H",
  y_comma = FALSE
)
save_panel(out$plot, file.path(RPT, "panels", "panel_h_1rm_ext"), 150, 100)
write_csv(out$audit, file.path(DAT, "audit_panel_H_1rm_ext.csv"))
