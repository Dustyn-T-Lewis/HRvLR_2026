if (!exists("meta")) source(here::here("04_Figures", "F01", "a_script", "setup.R"))

c_continuum <- build_continuum(meta, tag = "C")
save_panel(c_continuum$plot, file.path(F01_RPT, "panels", "panel_c_continuum"), 150, 120)
F01_PANELS[["continuum"]] <- c_continuum$plot
F01_AUDIT[["responder_axis"]] <- c_continuum$audit
