# F03 ring-volcano grid (enrichVolcano): one volcano-in-ring per contrast, all 9
# in a 3x3 grid. Volcano = DEP logFC/P; ring arcs = fgsea NES (top 7 per
# direction) with thickness = -log10(padj). Cache + DEP read from setup.
pacman::p_load(here, dplyr, readr, tibble, enrichVolcano, ggplot2, patchwork, openxlsx)
if (!exists("dep")) source(here("04_Figures", "F03", "a_script", "HRvLR_F03_setup.R"))

# Named lists (one frame per contrast) so a contrast with no significant
# enrichment (e.g. Trained_HRvLR) still draws its volcano with a bare ring.
volc_dfs <- lapply(CONTRASTS, function(ct) {
  tibble(
    gene = dep$gene,
    logFC = dep[[paste0("logFC_", ct)]],
    P.Value = dep[[paste0("P.Value_", ct)]],
    padj = dep[[paste0("adj.P.Val_", ct)]]
  ) |>
    filter(!is.na(logFC) & !is.na(P.Value))
})
names(volc_dfs) <- CONTRASTS

# Enrichment: significant terms, top 7 per direction, names cleaned (0-row ok)
enrich_dfs <- lapply(CONTRASTS, function(ct) {
  base <- fg |> filter(contrast == ct, is.finite(NES))
  e <- base |>
    filter(!is.na(padj), padj < 0.05) |>
    mutate(dir = if_else(NES > 0, "up", "dn")) |>
    group_by(dir) |>
    slice_min(padj, n = 7, with_ties = FALSE) |>
    ungroup()
  if (nrow(e) == 0) e <- slice_min(base, padj, n = 1) # no sig term: show the best one
  tibble(
    pathway = clean_pathway_name(e$pathway),
    NES = e$NES, padj = e$padj, size = e$size,
    leading_edge = e$leadingEdge, database = e$database
  )
})
names(enrich_dfs) <- CONTRASTS

g <- suppressMessages(volcano_ring_grid(
  volc_dfs, enrich_dfs,
  contrasts = CONTRASTS, ncol = 3, genes_sep = ";"
))

# The 3x3 ring-volcano grid is the whole figure, so it is the composite.
ggsave(file.path(RPT_DIR, "F03_composite.png"), g$plot,
  width = 380, height = 380, units = "mm", dpi = 200, bg = "white"
)
ggsave(file.path(RPT_DIR, "F03_composite.pdf"), g$plot,
  width = 380, height = 380, units = "mm", device = PDF_DEVICE, bg = "white"
)

# Source-data workbook: the significant enrichment behind every ring
sig <- fg |>
  filter(!is.na(padj), padj < 0.05, is.finite(NES)) |>
  transmute(
    contrast, database,
    pathway = clean_pathway_name(pathway),
    NES, padj, size, leading_edge = leadingEdge
  ) |>
  arrange(contrast, padj)
wb <- createWorkbook()
addWorksheet(wb, "fgsea_significant")
writeData(wb, "fgsea_significant", sig)
addWorksheet(wb, "fgsea_all")
writeData(wb, "fgsea_all", fg)
saveWorkbook(wb, file.path(DAT_DIR, "F03_source_data.xlsx"), overwrite = TRUE)
cat("F03 composite + workbook done.\n")
