# F03 ring inputs. Read the cached fgsea (fg), keep each contrast's significant
# pathways, collapse them with the EnrichmentMap display dedup, take the top
# RING_N by FDR, and shape the volcano and enrichment frames enrichVolcano draws.
# The cache (fgsea_all) is never touched here; all reduction is display-side.
pacman::p_load(dplyr, tibble, enrichVolcano, ggplot2, patchwork, shadowtext)

RING_N <- 12

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

# The `padj` slot is enrichVolcano's significance channel; we feed it the Xiao
# Eq.2 transformed P-value Pi = p^|log2FC| (not an adjusted p, not FDR). Display
# it as Pi, never as adjusted p.
ring_volc <- function(dep, ct) {
  tibble(
    gene = dep$gene,
    logFC = dep[[paste0("logFC_", ct)]],
    P.Value = dep[[paste0("P.Value_", ct)]],
    padj = dep[[paste0("pi_score_", ct)]]
  ) |>
    filter(!is.na(logFC), !is.na(P.Value))
}

# Flag the RING_N most-significant kept pathways as drawn; the report is already
# padj-ordered, so the first RING_N kept rows are the arcs on the ring.
mark_drawn <- function(report, ring_n = RING_N) {
  kept <- which(report$dedup_status == "kept")
  report$drawn <- FALSE
  report$drawn[kept[seq_len(min(ring_n, length(kept)))]] <- TRUE
  report
}

# One contrast as a single volcano_ring, tagged and recoloured to its family
# palette. Volcano points gate at PI_THRESH (Xiao Eq.2 Pi); ring arcs carry the
# fgsea BH q. Returns the plot and its dedup report for the workbook.
ring_plot <- function(fg, dep, pw, ct, palette, tag = NULL) {
  prep <- ring_enrich(fg, ct, pw)
  p <- suppressMessages(volcano_ring(
    ring_volc(dep, ct), prep$enrich,
    title = ct, tag = tag, genes_sep = ";", p_threshold = PI_THRESH,
    theme = volcano_ring_theme(
      up = palette$up, down = palette$down, nes_colors = palette$nes
    )
  ))
  list(plot = p, report = mark_drawn(prep$report))
}

# Assemble tagged rings into one composite. `specs` is an ordered list of
# list(contrast, palette, tag); `byrow = FALSE` fills down columns so the tags
# read A,C,E across the top and B,D,F below. guides = "collect" merges each
# family's identical NES scale into one legend.
build_ring_grid <- function(fg, dep, pw, specs, ncol, byrow = TRUE) {
  panels <- lapply(specs, function(s) {
    ring_plot(fg, dep, pw, s$contrast, s$palette, s$tag)
  })
  grid <- patchwork::wrap_plots(
    lapply(panels, `[[`, "plot"),
    ncol = ncol, byrow = byrow, guides = "collect"
  ) & ggplot2::theme(
    legend.position = "bottom",
    plot.margin = ggplot2::margin(2, 6, 2, 6, "mm")
  )
  list(plot = grid, reports = bind_rows(lapply(panels, `[[`, "report")))
}
