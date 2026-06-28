# Stitch the F04 composite from the rendered panels (A concordance map + flanking
# ORA, B pathway NES map) and write the single source-data workbook.
pacman::p_load(patchwork, openxlsx)

composite <- (wrap_elements(full = panel_a$plot) / panel_b$plot) +
  plot_layout(heights = c(1.5, 1.1)) +
  plot_annotation(
    title = "Figure 4. Training response: shared adaptation and responder-specific divergence",
    subtitle = "HR vs LR over T1 to T2. Diagonal = concordant; off-diagonal and orange = interaction-significant divergence.",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 2),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
    )
  )

save_panel(composite, file.path(RPT_DIR, "F04_composite"), 340, 320)

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
    "Panel A: per-protein HR/LR training logFC, quadrant, divergence class",
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
saveWorkbook(wb, file.path(DAT_DIR, "F04_source_data.xlsx"), overwrite = TRUE)
cat("F04 composite + workbook written.\n")
