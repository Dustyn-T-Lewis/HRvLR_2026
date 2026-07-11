if (!exists("meta")) source(here::here("04_Figures", "F01", "a_script", "setup.R"))

b_composite <- build_composite(meta, tag = "B")
save_panel(b_composite$plot, file.path(F01_RPT, "panels", "panel_b_composite"), 150, 100)
F01_PANELS[["composite"]] <- b_composite$plot
F01_AUDIT[["composite"]] <- b_composite$audit
