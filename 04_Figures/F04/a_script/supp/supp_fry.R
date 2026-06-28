# Supp - fry running-ES: do the HR-significant DEP sets move coordinately in LR?
if (!exists("dep")) {
  source(here::here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))
}

fry_out <- run_fry_concordance(da, dep, pw, C_HI, C_LO, LO_LEVELS)
supp_fry_plot <- panel_fry(fry_out, LABELS)

save_panel(supp_fry_plot, file.path(RPT_DIR, "supp", "supp_fry"), 185, 150)
write_csv(fry_out$results, file.path(DAT_DIR, "supp_fry_results.csv"))
