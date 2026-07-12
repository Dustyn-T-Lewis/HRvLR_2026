if (!exists("meta")) source(here::here("04_Figures", "F01", "a_script", "setup.R"))

c_forest <- build_forest(meta, tag = "C")
save_panel(c_forest$plot, file.path(F01_RPT, "panels", "panel_c_forest"), 150, 120)
F01_PANELS[["forest"]] <- c_forest$plot
F01_AUDIT[["change_advantage"]] <- c_forest$audit
