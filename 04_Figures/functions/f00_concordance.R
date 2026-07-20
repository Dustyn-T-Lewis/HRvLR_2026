# Driver for the HR-vs-LR concordance figures (F03_pathway/supp training and acute). One
# parameterised render keeps the two figures a single source of truth; only the
# per-figure config (contrast pair, factor levels, labels, composite titles)
# differs. Builds the 5-panel YvO engine composite: A quadrant ORA, B pattern
# heatmap + Sankey, C fry, D pathway NES scatter, E RRHO2. Reads the DEP results,
# the imputed matrix (fry), and the shared F03 fgsea cache (NES scatter, computed
# once in F03, never here).
pacman::p_load(here, dplyr, readr, tibble, patchwork, cowplot, ggplot2, openxlsx)

render_concordance_figure <- function(cfg) {
  source(here("functions", "shared_style.R"))
  source(here("functions", "shared_pathway_utils.R"))
  source(here("04_Figures", "functions", "f00_concordance_panels.R"))

  fig_dir <- if (is.null(cfg$fig_dir)) here("04_Figures", cfg$fig_id) else cfg$fig_dir
  rpt <- file.path(fig_dir, "b_reports")
  dat <- file.path(fig_dir, "c_data")
  panels <- file.path(rpt, "panels")
  supp <- file.path(rpt, "supp")
  for (d in c(panels, supp, dat)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  dep <- read_csv(
    here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
    show_col_types = FALSE
  )
  da <- readRDS(here("02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"))
  cache <- as_tibble(openxlsx::read.xlsx(
    here("04_Figures", "F03_pathway", "c_data", "F03_pathway_source_data.xlsx"),
    sheet = "fgsea_all"
  ))
  pw <- build_pathway_collection(
    min_size = 15, max_size = 500, include_goslim = FALSE, exclude_variants = TRUE
  )

  message("=== ", cfg$fig_id, " Panel A (quadrant ORA) ===")
  quad_tbl <- build_quadrant_table(dep, cfg$c_hi, cfg$c_lo, cfg$c_int)
  quad_ora <- run_quadrant_ora(quad_tbl, pw)
  pa <- panel_quadrant_ora(quad_tbl, quad_ora, cfg)
  save_panel(pa$plot, file.path(panels, "panel_a_quadrant_ora"), 250, 165)

  message("=== ", cfg$fig_id, " Panel B (pattern heatmap) ===")
  pb <- panel_pattern_heatmap(dep, cfg)
  save_panel(pb$plot, file.path(panels, "panel_b_pattern_heatmap"), 200, 150)

  message("=== ", cfg$fig_id, " Panel C (fry) ===")
  fry_out <- run_fry_concordance(da, dep, pw, cfg$c_hi, cfg$c_lo, cfg$lo_levels)
  pc <- panel_fry(fry_out, cfg)
  save_panel(pc, file.path(panels, "panel_c_fry"), 220, 170)

  message("=== ", cfg$fig_id, " Panel D (NES scatter) ===")
  nes_wide <- build_nes_wide(cache, cfg$c_hi, cfg$c_lo, cfg$c_int)
  pd <- panel_nes_scatter(nes_wide, cfg$c_hi, cfg$c_lo, cfg)
  save_panel(pd$plot, file.path(panels, "panel_d_nes_scatter"), 150, 160)

  message("=== ", cfg$fig_id, " Panel E (RRHO2) ===")
  pe <- panel_rrho2(dep, cfg$c_hi, cfg$c_lo, cfg)
  save_panel(pe$plot, file.path(panels, "panel_e_rrho2"), 120, 120)

  # Supplements
  message("=== ", cfg$fig_id, " supplements ===")
  supp_thresh <- panel_threshold_sens(dep, cfg$c_hi, cfg$c_lo)
  save_panel(supp_thresh$plot, file.path(supp, "supp_threshold"), 110, 90)
  supp_boot <- panel_rho_bootstrap(dep, cfg$c_hi, cfg$c_lo)
  save_panel(supp_boot$plot, file.path(supp, "supp_bootstrap"), 110, 90)
  fry_le <- (ora_bars(fry_out$ora_up, COMP_RED) / ora_bars(fry_out$ora_down, COMP_BLUE)) +
    plot_annotation(
      title = "fry leading-edge ORA",
      subtitle = "Top concordant driving proteins per direction",
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
      )
    )
  save_panel(fry_le, file.path(supp, "supp_fry_leading"), 150, 130)

  # Composite: A + B top, C + D + E bottom (YvO F04 geometry)
  message("=== ", cfg$fig_id, " composite ===")
  ttl_A <- "Quadrant ORA (Concordance)"
  sub_A <- sprintf(
    "N = %d | %d DEPs (\u03a0) | %d enriched (FDR) | \u03c1 = %.2f",
    pa$n_total, pa$n_sig, pa$n_enrich, pa$rho
  )
  ttl_B <- "Protein-to-Pathway"
  sub_B <- sprintf("%d proteins | %d GO-Slim categories", pb$n_total, pb$n_pw)
  ttl_C <- "fry: Concordance"
  sub_C <- sprintf("n = %d proteins | dupCor = %.3f", fry_out$n_all, fry_out$cor_within)
  ttl_D <- "Pathway Concordance"
  sub_D <- sprintf("\u03c1 = %.2f | %.0f%% concordant", pd$rho, pd$conc_frac * 100)
  ttl_E <- "RRHO2 Concordance"
  sub_E <- sprintf("%d genes | %d concordant hotspot", pe$n_shared, pe$n_concordant)

  pd_comp <- pd$plot + labs(subtitle = NULL)
  pe_comp <- pe$plot

  layout <- paste(
    "##############",
    "AAAAAAAABBBBBB",
    "AAAAAAAABBBBBB",
    "AAAAAAAABBBBBB",
    "AAAAAAAABBBBBB",
    "AAAAAAAABBBBBB",
    "AAAAAAAABBBBBB",
    "##############",
    "##############",
    "CCCCCCDDDDEEEE",
    "CCCCCCDDDDEEEE",
    "CCCCCCDDDDEEEE",
    "CCCCCCDDDDEEEE",
    "CCCCCCDDDDEEEE",
    "CCCCCCDDDDEEEE",
    sep = "\n"
  )

  fig <- wrap_elements(full = pa$plot) +
    wrap_elements(full = pb$plot) +
    wrap_elements(full = pc) +
    wrap_elements(full = pd_comp) +
    wrap_elements(full = pe_comp) +
    plot_layout(
      design = layout,
      widths = rep(1, 14),
      heights = c(6.5, rep(10, 6), 4, 4.5, rep(12, 6))
    )

  TAG_SZ <- 20
  TTL_SZ <- 15
  SUB_SZ <- 10
  X_TTL <- 0.028
  SUB_OFF <- 0.020

  draw_head <- function(g, tag, ttl, sub, x, y) {
    g +
      draw_label(tag, x = x, y = y, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
      draw_label(ttl, x = x + X_TTL, y = y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
      draw_label(sub,
        x = x + X_TTL, y = y - SUB_OFF, size = SUB_SZ, fontface = "bold.italic",
        hjust = 0, vjust = 1, colour = "grey40"
      )
  }

  composite <- ggdraw(fig)
  composite <- draw_head(composite, "A", ttl_A, sub_A, 0.006, 0.992)
  composite <- draw_head(composite, "B", ttl_B, sub_B, 0.560, 0.992)
  composite <- draw_head(composite, "C", ttl_C, sub_C, 0.006, 0.520)
  composite <- draw_head(composite, "D", ttl_D, sub_D, 0.430, 0.520)
  composite <- draw_head(composite, "E", ttl_E, sub_E, 0.720, 0.520)
  composite <- composite +
    draw_label(
      "Panel D (NES concordance) uses Hallmark + GO Slim; ORA panels use the full pathway collection.",
      x = 0.5, y = 0.006, size = SUB_SZ * 0.85, fontface = "italic",
      hjust = 0.5, vjust = 0, colour = "grey45"
    )

  ggsave(file.path(rpt, paste0(cfg$fig_id, "_composite.pdf")), composite,
    width = 420, height = 310, units = "mm", device = PDF_DEVICE
  )
  ggsave(file.path(rpt, paste0(cfg$fig_id, "_composite.png")), composite,
    width = 420, height = 310, units = "mm", dpi = 300
  )

  # Source-data workbook
  sheets <- list(
    panel_a_quadrants = quad_tbl,
    panel_a_ora = quad_ora,
    panel_b_pattern = pb$data |> select(gene, quadrant, sig_cat, pathway, lfc_x, lfc_y),
    panel_c_fry_results = fry_out$results,
    panel_d_nes = pd$data,
    panel_e_rrho2_hotspot = pe$hotspot,
    supp_threshold = supp_thresh$data,
    supp_bootstrap = supp_boot$data
  )
  overview <- data.frame(
    sheet = names(sheets),
    description = c(
      sprintf("Panel A: per-protein HR/LR %s logFC, quadrant, sig class", tolower(cfg$labels$phase)),
      "Panel A: over-representation per concordance quadrant",
      "Panel B: pi-sig proteins, quadrant, GO-Slim category, HR/LR logFC",
      "Panel C: fry rotation test of HR DEP sets along the LR contrast",
      "Panel D: pathway NES (HR vs LR) with cached interaction NES/padj",
      "Panel E: RRHO2 hotspot genes per corner",
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

  message(sprintf("%s complete: %s vs %s", cfg$fig_id, cfg$c_hi, cfg$c_lo))
}
