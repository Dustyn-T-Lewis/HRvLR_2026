if (!exists("meta")) source(here::here("04_Figures", "F01", "a_script", "setup.R"))

b_continuum <- build_continuum(meta, tag = "B")
save_panel(b_continuum$plot, file.path(F01_RPT, "panels", "panel_b_continuum"), 150, 120)
F01_PANELS[["continuum"]] <- b_continuum$plot
F01_AUDIT[["responder_axis"]] <- b_continuum$audit
