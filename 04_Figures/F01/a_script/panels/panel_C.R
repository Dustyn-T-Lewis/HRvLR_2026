# F01 Panel C: fCSA Mixed (Pre/Post x Group + Delta)
pacman::p_load(here, readr)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "F01", "a_script", "functions", "prepost.R"))

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

out <- prepost_delta_panel(
  "fCSA_Mixed_Pre", "fCSA Mixed",
  expression("fCSA (" * mu * m^2 * ")"), expression(Delta ~ "fCSA (" * mu * m^2 * ")"),
  "C"
)
save_panel(out$plot, file.path(RPT, "panels", "panel_c_fcsa_mixed"), 150, 100)
write_csv(out$audit, file.path(DAT, "audit_panel_C_fcsa_mixed.csv"))
