# Supp - RRHO2 threshold-free overlap of the HR and LR training ranks.
if (!exists("dep")) {
  source(here::here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))
}

rrho <- panel_rrho2(dep, C_HI, C_LO, LABELS)

save_panel(rrho$plot, file.path(RPT_DIR, "supp", "supp_rrho2"), 110, 110)
