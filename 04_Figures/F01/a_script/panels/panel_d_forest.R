if (!exists("meta")) source(here::here("04_Figures", "F01", "a_script", "setup.R"))

d_forest <- build_forest(meta, tag = "D")
save_panel(d_forest$plot, file.path(F01_RPT, "panels", "panel_d_forest"), 150, 120)
F01_PANELS[["forest"]] <- d_forest$plot
F01_AUDIT[["lmm_interaction"]] <- d_forest$audit
