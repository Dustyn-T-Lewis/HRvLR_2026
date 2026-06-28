# Supp - quadrant split stability across significance thresholds.
if (!exists("dep")) {
  source(here::here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))
}

supp_thresh <- panel_threshold_sens(dep, C_HI, C_LO)

save_panel(supp_thresh$plot, file.path(RPT_DIR, "supp", "supp_threshold"), 110, 90)
write_csv(supp_thresh$data, file.path(DAT_DIR, "supp_threshold.csv"))
