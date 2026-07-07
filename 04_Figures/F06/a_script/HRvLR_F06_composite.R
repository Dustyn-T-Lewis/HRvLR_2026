# F06 composite: stitch the rendered module panels (A atlas, B ORA, C phenotype)
# and write the source-data workbook. Run after the panels exist.
pacman::p_load(here, patchwork, ggplot2, png, grid, readr, openxlsx)
source(here("04_Figures", "functions", "style.R"))

RPT_DIR <- here("04_Figures", "F06", "b_reports")
DAT_DIR <- here("04_Figures", "F06", "c_data")

read_panel <- function(name) {
  path <- file.path(RPT_DIR, "panels", paste0(name, ".png"))
  if (!file.exists(path)) stop("Missing panel: ", path, " (run HRvLR_F06_run.R first)")
  ggplot() +
    annotation_custom(rasterGrob(readPNG(path), interpolate = TRUE)) +
    theme_void()
}

composite <- (read_panel("panel_a_wgcna_map") /
  (read_panel("panel_b_wgcna_ora") | read_panel("panel_c_module_phenotype"))) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title = "Figure 6. Co-expression module atlas and phenotype links",
    subtitle = "Unsupervised WGCNA modules, their enriched biology, and how each module's response tracks the phenotype.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT_DIR, "F06_composite.png"), composite,
  width = 460, height = 420, units = "mm", dpi = 200, bg = "white"
)
ggsave(file.path(RPT_DIR, "F06_composite.pdf"), composite,
  width = 460, height = 420, units = "mm", device = PDF_DEVICE, bg = "white"
)

sheet_desc <- c(
  wgcna_membership = "Protein-to-module assignment (WGCNA, full imputed proteome)",
  wgcna_eigengene = "Module eigengenes per sample and timepoint",
  module_prediction = "Baseline (T1) eigengene LOO-CV prediction of adaptation traits",
  module_responder = "Per-module responder (HR vs LR) signal per timepoint",
  module_ora = "Panel B: over-representation per module",
  module_phenotype = "Panel C: module eigengene delta vs response trait",
  module_phenotype_lmm = "Panel S: module-phenotype linear mixed model",
  phenotype = "Per-subject phenotype table (T2 composite + T1->T2 deltas)"
)
present <- names(sheet_desc)[file.exists(
  file.path(DAT_DIR, paste0(names(sheet_desc), ".csv"))
)]
overview <- data.frame(
  sheet = substr(present, 1, 31), description = unname(sheet_desc[present]),
  stringsAsFactors = FALSE
)
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
for (s in present) {
  addWorksheet(wb, substr(s, 1, 31))
  writeData(wb, substr(s, 1, 31), read_csv(
    file.path(DAT_DIR, paste0(s, ".csv")),
    show_col_types = FALSE
  ))
}
saveWorkbook(wb, file.path(DAT_DIR, "F06_source_data.xlsx"), overwrite = TRUE)
cat("F06 composite + workbook saved to", RPT_DIR, "\n")
