# Supp - bootstrap CI on the HR-vs-LR logFC Spearman concordance.
if (!exists("dep")) {
  source(here::here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))
}

supp_boot <- panel_rho_bootstrap(dep, C_HI, C_LO)

save_panel(supp_boot$plot, file.path(RPT_DIR, "supp", "supp_bootstrap"), 110, 90)
write_csv(supp_boot$data, file.path(DAT_DIR, "supp_bootstrap.csv"))
