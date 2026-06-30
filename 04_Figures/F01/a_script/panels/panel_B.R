# F01 Panel B: 1RM Leg Press (Pre/Post x Group + Delta)
pacman::p_load(here, readr)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "F01", "a_script", "functions", "prepost.R"))

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

out <- prepost_delta_panel(
  "X1RM_Leg_Pre", "1RM Leg Press", "1RM (kg)", "Δ 1RM (kg)", "B",
  y_comma = FALSE
)
save_panel(out$plot, file.path(RPT, "panels", "panel_b_1rm_leg"), 150, 100)
write_csv(out$audit, file.path(DAT, "audit_panel_B_1rm_leg.csv"))
