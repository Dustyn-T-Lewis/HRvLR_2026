# F03 ring-volcanoes: one volcano-in-ring per contrast (volcano = DEP logFC vs
# -log10 P, points coloured by the pipeline pi-gate pi_score < 0.05; ring arcs =
# fgsea NES, top 7 per direction, thickness = -log10 padj). Two composites: the
# within-group HR/LR responses, and the between-responder + interaction contrasts.
pacman::p_load(here, dplyr, readr, tibble, enrichVolcano, ggplot2, patchwork)
if (!exists("dep")) source(here("04_Figures", "F03", "a_script", "HRvLR_F03_setup.R"))

volc_dfs <- lapply(CONTRASTS, function(ct) {
  tibble(
    gene = dep$gene,
    logFC = dep[[paste0("logFC_", ct)]],
    P.Value = dep[[paste0("P.Value_", ct)]],
    padj = dep[[paste0("pi_score_", ct)]]
  ) |>
    filter(!is.na(logFC) & !is.na(P.Value))
})
names(volc_dfs) <- CONTRASTS

enrich_dfs <- lapply(CONTRASTS, function(ct) {
  base <- fg |> filter(contrast == ct, is.finite(NES))
  e <- base |>
    filter(!is.na(padj), padj < 0.05) |>
    mutate(dir = if_else(NES > 0, "up", "dn")) |>
    group_by(dir) |>
    slice_min(padj, n = 7, with_ties = FALSE) |>
    ungroup()
  if (nrow(e) == 0) e <- slice_min(base, padj, n = 1)
  tibble(
    pathway = clean_pathway_name(e$pathway),
    NES = e$NES, padj = e$padj, size = e$size,
    leading_edge = e$leadingEdge, database = e$database
  )
})
names(enrich_dfs) <- CONTRASTS

save_grid <- function(contrasts, stem, ncol, width, height) {
  g <- suppressMessages(volcano_ring_grid(
    volc_dfs, enrich_dfs,
    contrasts = contrasts, ncol = ncol, genes_sep = ";"
  ))
  ggsave(file.path(RPT_DIR, paste0(stem, ".png")), g$plot,
    width = width, height = height, units = "mm", dpi = 200, bg = "white"
  )
  ggsave(file.path(RPT_DIR, paste0(stem, ".pdf")), g$plot,
    width = width, height = height, units = "mm", device = PDF_DEVICE, bg = "white"
  )
}

save_grid(c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR"),
  "F03_responses",
  ncol = 2, width = 380, height = 380
)
save_grid(
  c(
    "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
    "Training_Interaction", "Acute_Interaction"
  ),
  "F03_differential",
  ncol = 3, width = 570, height = 420
)

cat("F03 two composites done.\n")
