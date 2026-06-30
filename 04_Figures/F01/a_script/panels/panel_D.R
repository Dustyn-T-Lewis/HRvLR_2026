# F01 Panel D: fCSA Type I (Pre/Post x Group + Delta)
pacman::p_load(here, readr)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "F01", "a_script", "functions", "prepost.R"))

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

out <- prepost_delta_panel(
  "fCSA_Type_I_Pre", "fCSA Type I",
  expression("fCSA (" * mu * m^2 * ")"), expression(Delta ~ "fCSA (" * mu * m^2 * ")"),
  "D"
)
save_panel(out$plot, file.path(RPT, "panels", "panel_d_fcsa_type1"), 150, 100)
write_csv(out$audit, file.path(DAT, "audit_panel_D_fcsa_type1.csv"))
