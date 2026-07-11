# F03 ring inputs. Read the cached fgsea (fg), keep each contrast's significant
# pathways, collapse them with the EnrichmentMap display dedup, take the top
# RING_N by FDR, and shape the volcano and enrichment frames enrichVolcano draws.
# The cache (fgsea_all) is never touched here; all reduction is display-side.
pacman::p_load(dplyr, tibble, enrichVolcano, ggplot2, patchwork)

RING_N <- 12
TOP30_N <- 30

ring_significant <- function(fg, ct) {
  fg |>
    filter(contrast == ct, !is.na(padj), padj < 0.05, is.finite(NES)) |>
    arrange(padj)
}

# Significant -> EnrichmentMap dedup (report) -> top RING_N kept by FDR. Falls
# back to the single best pathway when a contrast has no significant enrichment,
# so the ring still draws.
ring_enrich <- function(fg, ct, pw, ring_n = RING_N) {
  sig <- as.data.frame(ring_significant(fg, ct))
  if (nrow(sig) == 0) {
    sig <- fg |>
      filter(contrast == ct, is.finite(NES)) |>
      slice_min(padj, n = 1) |>
      as.data.frame()
  }
  report <- dedup_report(sig, pw)
  kept <- report[report$dedup_status == "kept", , drop = FALSE]
  kept <- kept[order(kept$padj), , drop = FALSE]
  kept <- kept[seq_len(min(ring_n, nrow(kept))), , drop = FALSE]
  enrich <- tibble(
    pathway = clean_pathway_name(kept$pathway),
    NES = kept$NES, padj = kept$padj, size = kept$size,
    leading_edge = kept$leadingEdge, database = kept$database
  )
  list(enrich = enrich, report = as_tibble(report))
}

ring_volc <- function(dep, ct) {
  tibble(
    gene = dep$gene,
    logFC = dep[[paste0("logFC_", ct)]],
    P.Value = dep[[paste0("P.Value_", ct)]],
    padj = dep[[paste0("pi_score_", ct)]]
  ) |>
    filter(!is.na(logFC), !is.na(P.Value))
}

# No-dedup top-N up and top-N down by FDR, for the supplement bars.
top30_updown <- function(fg, ct, n = TOP30_N) {
  base <- fg |> filter(contrast == ct, !is.na(padj), is.finite(NES))
  up <- base |>
    filter(NES > 0) |>
    slice_min(padj, n = n, with_ties = FALSE)
  dn <- base |>
    filter(NES < 0) |>
    slice_min(padj, n = n, with_ties = FALSE)
  bind_rows(up, dn) |>
    transmute(
      contrast, database,
      pathway = clean_pathway_name(pathway),
      direction = if_else(NES > 0, "up", "down"),
      NES, padj, size
    )
}

# Flag the RING_N most-significant kept pathways as drawn; the report is already
# padj-ordered, so the first RING_N kept rows are the arcs on the ring.
mark_drawn <- function(report, ring_n = RING_N) {
  kept <- which(report$dedup_status == "kept")
  report$drawn <- FALSE
  report$drawn[kept[seq_len(min(ring_n, length(kept)))]] <- TRUE
  report
}

# One contrast family as a volcano_ring_grid, recoloured to the family palette.
# Returns the grid object and the per-contrast dedup reports for the workbook.
render_family <- function(fg, dep, pw, contrasts, palette, ncol) {
  prep <- lapply(contrasts, function(ct) ring_enrich(fg, ct, pw))
  names(prep) <- contrasts
  volc_dfs <- lapply(contrasts, function(ct) ring_volc(dep, ct))
  names(volc_dfs) <- contrasts
  enrich_dfs <- lapply(prep, `[[`, "enrich")
  grid <- suppressMessages(volcano_ring_grid(
    volc_dfs, enrich_dfs,
    contrasts = contrasts, ncol = ncol, genes_sep = ";",
    tag_levels = NULL, legend_position = "bottom",
    theme = volcano_ring_theme(
      up = palette$up, down = palette$down, nes_colors = palette$nes
    )
  ))
  reports <- Map(function(r, ct) mark_drawn(r[["report"]]), prep, contrasts)
  list(grid = grid, reports = bind_rows(reports))
}
