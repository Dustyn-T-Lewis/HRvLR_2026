# Supplement: the no-dedup top-30 up and top-30 down pathways by FDR for each of
# the seven main-composite contrasts, so the rings' parsimony is auditable against
# the full ranked enrichment. Bars point up (right) or down (left) and are filled
# by database; each pathway name runs along its bar, freeing the label gutter.
MAIN_CONTRASTS <- c(
  "Training_HR", "Training_LR", "Acute_HR", "Acute_LR",
  "Baseline_HRvLR", "Training_Interaction", "Acute_Interaction"
)

top30_bar <- function(tab, ct) {
  d <- tab |>
    filter(contrast == ct) |>
    mutate(
      pathway = stats::reorder(pathway, NES),
      label_x = if_else(NES > 0, 0.02, -0.02),
      label_h = if_else(NES > 0, 0, 1)
    )
  m <- max(abs(d$NES)) * 1.55
  ggplot(d, aes(NES, pathway, fill = database)) +
    geom_col(width = 0.85) +
    geom_text(aes(x = label_x, label = pathway, hjust = label_h),
      size = 2.4, colour = "grey10"
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
    scale_fill_manual(values = DB_COLORS, drop = FALSE) +
    scale_x_continuous(limits = c(-m, m)) +
    labs(title = ct, x = "NES", y = NULL, fill = "Database") +
    FIG_THEME +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

panel_top30 <- function(fg) {
  tab <- bind_rows(lapply(MAIN_CONTRASTS, function(ct) top30_updown(fg, ct)))
  plots <- lapply(MAIN_CONTRASTS, function(ct) top30_bar(tab, ct))
  combined <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  list(plot = combined, table = tab)
}
