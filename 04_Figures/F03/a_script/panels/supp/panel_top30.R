# Supplement: one figure per contrast, the top-30 down- and up-regulated pathways
# by NES sitting back to back and aligned strongest-at-top, no dedup, so the rings'
# parsimony is auditable against the full ranked enrichment. Bars are filled by
# database with the name inside in white bold and the FDR at each tip.
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
