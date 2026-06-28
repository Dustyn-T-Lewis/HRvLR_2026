# Panel B - pathway NES scatter (Hallmark + GO Slim), divergent pathways flagged.
if (!exists("dep")) {
  source(here::here("04_Figures", "F05", "a_script", "HRvLR_F05_setup.R"))
}

nes_wide <- build_nes_wide(cache, C_HI, C_LO, C_INT)
panel_b <- panel_nes_scatter(nes_wide, C_HI, C_LO, LABELS)

save_panel(panel_b$plot, file.path(RPT_DIR, "panels", "panel_b_nes_scatter"), 150, 150)
write_csv(panel_b$data, file.path(DAT_DIR, "panel_b_nes_scatter.csv"))
