# Supp - RRHO2 threshold-free overlap of the HR and LR acute ranks.
if (!exists("dep")) {
  source(here::here("04_Figures", "F05", "a_script", "HRvLR_F05_setup.R"))
}

rrho <- panel_rrho2(dep, C_HI, C_LO, LABELS)

save_panel(rrho$plot, file.path(RPT_DIR, "supp", "supp_rrho2"), 110, 110)
