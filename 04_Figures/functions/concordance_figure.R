# Driver for the HR-vs-LR concordance figures (F04 training, F05 acute). One
# parameterised render keeps the two figures a single source of truth; only the
# per-figure config (contrast pair, factor levels, labels, composite titles)
# differs. Reads the DEP results, the imputed matrix (fry), and the shared F03
# fgsea cache (NES scatter - computed once in F03, never here).
pacman::p_load(here, dplyr, readr, tibble, patchwork, openxlsx)

render_concordance_figure <- function(cfg) {
  source(here("04_Figures", "functions", "style.R"))
  source(here("04_Figures", "functions", "pathway_utils.R"))
  source(here("04_Figures", "functions", "concordance.R"))

  rpt <- here("04_Figures", cfg$fig_id, "b_reports")
  dat <- here("04_Figures", cfg$fig_id, "c_data")
  dir.create(file.path(rpt, "panels"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(rpt, "supp"), recursive = TRUE, showWarnings = FALSE)
  dir.create(dat, recursive = TRUE, showWarnings = FALSE)

  dep <- read_csv(
    here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
    show_col_types = FALSE
  )
  da <- readRDS(here("02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"))
  cache <- as_tibble(openxlsx::read.xlsx(
    here("04_Figures", "F03", "c_data", "F03_source_data.xlsx"),
    sheet = "fgsea_all"
  ))
  pw <- build_pathway_collection(
    min_size = 15, max_size = 500, include_goslim = FALSE, exclude_variants = TRUE
  )

  quad_tbl <- build_quadrant_table(dep, cfg$c_hi, cfg$c_lo, cfg$c_int)
  quad_ora <- run_quadrant_ora(quad_tbl, pw)
  panel_a <- panel_quadrant_ora(quad_tbl, quad_ora, cfg$labels)
  save_panel(panel_a$plot, file.path(rpt, "panels", "panel_a_quadrant_ora"), 250, 150)
  write_csv(quad_tbl, file.path(dat, "panel_a_quadrant_table.csv"))
  write_csv(quad_ora, file.path(dat, "panel_a_ora_quadrant.csv"))

  nes_wide <- build_nes_wide(cache, cfg$c_hi, cfg$c_lo, cfg$c_int)
  panel_b <- panel_nes_scatter(nes_wide, cfg$c_hi, cfg$c_lo, cfg$labels)
  save_panel(panel_b$plot, file.path(rpt, "panels", "panel_b_nes_scatter"), 150, 150)
  write_csv(panel_b$data, file.path(dat, "panel_b_nes_scatter.csv"))

  fry_out <- run_fry_concordance(da, dep, pw, cfg$c_hi, cfg$c_lo, cfg$lo_levels)
  save_panel(panel_fry(fry_out, cfg$labels), file.path(rpt, "supp", "supp_fry"), 185, 150)
  write_csv(fry_out$results, file.path(dat, "supp_fry_results.csv"))

  rrho <- panel_rrho2(dep, cfg$c_hi, cfg$c_lo, cfg$labels)
  save_panel(rrho$plot, file.path(rpt, "supp", "supp_rrho2"), 110, 110)

  supp_thresh <- panel_threshold_sens(dep, cfg$c_hi, cfg$c_lo)
  save_panel(supp_thresh$plot, file.path(rpt, "supp", "supp_threshold"), 110, 90)
  write_csv(supp_thresh$data, file.path(dat, "supp_threshold.csv"))

  supp_boot <- panel_rho_bootstrap(dep, cfg$c_hi, cfg$c_lo)
  save_panel(supp_boot$plot, file.path(rpt, "supp", "supp_bootstrap"), 110, 90)
  write_csv(supp_boot$data, file.path(dat, "supp_bootstrap.csv"))

  composite <- (wrap_elements(full = panel_a$plot) / panel_b$plot) +
    plot_layout(heights = c(1.5, 1.1)) +
    plot_annotation(
      title = cfg$title, subtitle = cfg$subtitle, tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 2),
        plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
      )
    )
  save_panel(composite, file.path(rpt, paste0(cfg$fig_id, "_composite")), 340, 320)

  sheets <- list(
    panel_a_quadrants = quad_tbl,
    panel_a_ora = quad_ora,
    panel_b_nes = panel_b$data,
    supp_fry_results = fry_out$results,
    supp_threshold = supp_thresh$data,
    supp_bootstrap = supp_boot$data
  )
  overview <- data.frame(
    sheet = names(sheets),
    description = c(
      sprintf("Panel A: per-protein HR/LR %s logFC, quadrant, divergence class", tolower(cfg$labels$phase)),
      "Panel A: over-representation per concordance quadrant",
      "Panel B: pathway NES (HR vs LR) with cached interaction NES/padj",
      "Supp: fry rotation test of HR DEP sets along the LR contrast",
      "Supp: quadrant counts across significance thresholds",
      "Supp: Spearman rho bootstrap (observed + 95% CI)"
    ),
    stringsAsFactors = FALSE
  )
  wb <- createWorkbook()
  addWorksheet(wb, "overview")
  writeData(wb, "overview", overview)
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, as.data.frame(sheets[[nm]]))
  }
  saveWorkbook(wb, file.path(dat, paste0(cfg$fig_id, "_source_data.xlsx")), overwrite = TRUE)
  cat(sprintf("%s complete: %s vs %s\n", cfg$fig_id, cfg$c_hi, cfg$c_lo))
}
