# Panel A - HR vs LR training logFC scatter + flanking per-quadrant ORA.
if (!exists("dep")) {
  source(here::here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))
}

quad_tbl <- build_quadrant_table(dep, C_HI, C_LO, C_INT)
quad_ora <- run_quadrant_ora(quad_tbl, pw)
panel_a <- panel_quadrant_ora(quad_tbl, quad_ora, LABELS)

save_panel(panel_a$plot, file.path(RPT_DIR, "panels", "panel_a_quadrant_ora"), 250, 150)
write_csv(quad_tbl, file.path(DAT_DIR, "panel_a_quadrant_table.csv"))
write_csv(quad_ora, file.path(DAT_DIR, "panel_a_ora_quadrant.csv"))
