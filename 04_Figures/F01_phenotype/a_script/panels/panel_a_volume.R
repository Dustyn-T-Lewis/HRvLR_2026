if (!exists("meta")) source(here::here("04_Figures", "F01", "a_script", "setup.R"))

a_volume <- build_volume(meta, tag = "A")
save_panel(a_volume$plot, file.path(F01_RPT, "panels", "panel_a_volume"), 150, 100)
F01_PANELS[["volume"]] <- a_volume$plot
F01_AUDIT[["volume_load"]] <- a_volume$audit
