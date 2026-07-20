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

# Top-30 audit: the un-deduplicated ranked enrichment behind each ring, so the
# rings' parsimony can be checked against the full list.
TOP30_N <- 30L

top_pathways <- function(g, direction) {
  g <- if (direction == "up") filter(g, NES > 0) else filter(g, NES < 0)
  g |>
    slice_max(if (direction == "up") NES else -NES, n = TOP30_N, with_ties = FALSE) |>
    arrange(if (direction == "up") NES else desc(NES))
}

top30_side <- function(g, direction) {
  up <- direction == "up"
  d <- g |>
    mutate(
      rid = factor(paste(database, pathway), levels = paste(database, pathway)),
      name = clean_pathway_name(pathway, max_chars = 32),
      fdr = if_else(padj < 0.001,
        formatC(padj, format = "e", digits = 0),
        formatC(padj, format = "f", digits = 3)
      )
    )
  ggplot(d, aes(NES, rid, fill = database)) +
    geom_col(width = 0.8, colour = "black", linewidth = 0.25) +
    shadowtext::geom_shadowtext(aes(label = name),
      x = if (up) 0.05 else -0.05, hjust = if (up) 0 else 1,
      size = 2.5, fontface = "bold", colour = "white", bg.colour = "grey20"
    ) +
    geom_text(aes(label = fdr),
      hjust = if (up) -0.18 else 1.18,
      size = 1.95, fontface = "bold", colour = "grey20"
    ) +
    scale_fill_manual(values = DB_COLORS, name = NULL, drop = FALSE) +
    scale_y_discrete(labels = NULL) +
    scale_x_continuous(
      limits = if (up) c(0, max(d$NES)) else c(min(d$NES), 0),
      expand = expansion(mult = if (up) c(0.02, 0.16) else c(0.16, 0.02))
    ) +
    coord_cartesian(clip = "off") +
    labs(title = if (up) "Up-regulated" else "Down-regulated", x = "NES", y = NULL) +
    FIG_THEME +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

top30_figure <- function(gsea, ct) {
  g <- filter(gsea, contrast == ct)
  (top30_side(top_pathways(g, "down"), "down") |
    top30_side(top_pathways(g, "up"), "up")) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%s: top %d up / %d down by NES, no dedup", ct, TOP30_N, TOP30_N),
      theme = theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE))
    ) &
    theme(legend.position = "bottom")
}

panel_top30 <- function(fg, contrasts = CONTRASTS) {
  gsea <- fg |> filter(contrast %in% contrasts, is.finite(NES), !is.na(padj))
  plots <- lapply(contrasts, function(ct) top30_figure(gsea, ct))
  names(plots) <- contrasts
  table <- bind_rows(lapply(contrasts, function(ct) {
    g <- filter(gsea, contrast == ct)
    bind_rows(up = top_pathways(g, "up"), down = top_pathways(g, "down"), .id = "direction") |>
      transmute(
        contrast = ct, direction, database,
        pathway = clean_pathway_name(pathway),
        NES = round(NES, 2), padj = signif(padj, 3), size
      )
  }))
  list(plots = plots, table = table)
}
