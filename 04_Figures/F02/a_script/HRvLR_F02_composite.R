# F02 composite: six title-less main panels tiled 2x3 in a full-width YvO-style grid,
# with tags, titles, and subtitles drawn on the composite (not baked into the panels).
# Row 1: A state, B trajectory, C divergence; Row 2: D effect size, E DEPs, F pathways.
# Run HRvLR_F02_run.R first.
pacman::p_load(here, ggplot2, cowplot, png, grid)

source(here("04_Figures", "functions", "style.R"))

RPT_DIR <- here("04_Figures", "F02", "b_reports")
PANEL_DIR <- file.path(RPT_DIR, "panels")

read_panel <- function(stem) {
  path <- file.path(PANEL_DIR, paste0(stem, ".png"))
  if (!file.exists(path)) stop("Missing panel: ", path, " (run HRvLR_F02_run.R first)")
  rasterGrob(readPNG(path), interpolate = TRUE)
}

panels <- list(
  list(
    tag = "A", stem = "panel_a_pca", aspect = 140 / 95, row = 1, col = 1,
    title = "Global Proteome State", subtitle = "PCA — groups overlap (PERMANOVA p = 0.62)"
  ),
  list(
    tag = "B", stem = "panel_b_trajectory", aspect = 120 / 95, row = 1, col = 2,
    title = "Response Magnitude ‖Δ‖", subtitle = "RRPP — HR and LR move equally far"
  ),
  list(
    tag = "C", stem = "panel_c_divergence", aspect = 150 / 95, row = 1, col = 3,
    title = "Responder Divergence", subtitle = "HR vs LR logFC; interaction-significant labelled"
  ),
  list(
    tag = "D", stem = "panel_d_effectsize", aspect = 120 / 95, row = 2, col = 1,
    title = "Effect-Size Distribution", subtitle = "logFC per contrast"
  ),
  list(
    tag = "E", stem = "panel_e_dep_counts", aspect = 130 / 95, row = 2, col = 2,
    title = "DEPs per Contrast", subtitle = "% of proteome; dotted = chance"
  ),
  list(
    tag = "F", stem = "panel_f_pathways", aspect = 150 / 95, row = 2, col = 3,
    title = "Pathway Enrichment (GO)", subtitle = "non-redundant sets, BH < 0.05"
  )
)

MAR <- 6
GUT_X <- 5
GUT_Y <- 6
TITLE_H <- 9
COMP_W <- 185
col_w <- (COMP_W - 2 * MAR - 2 * GUT_X) / 3

img_h <- vapply(panels, function(p) col_w / p$aspect, numeric(1))
row_h <- c(
  max(img_h[vapply(panels, `[[`, numeric(1), "row") == 1]),
  max(img_h[vapply(panels, `[[`, numeric(1), "row") == 2])
)
cell_h <- TITLE_H + row_h
COMP_H <- sum(cell_h) + GUT_Y + 2 * MAR

col_x <- function(col) MAR + (col - 1) * (col_w + GUT_X)
title_top <- function(row) COMP_H - MAR - sum(cell_h[seq_len(row - 1)]) - GUT_Y * (row - 1)

TAG_SZ <- 8
TTL_SZ <- 7
SUB_SZ <- 5

composite <- ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA))
for (i in seq_along(panels)) {
  p <- panels[[i]]
  tx <- col_x(p$col)
  ty <- title_top(p$row)
  composite <- composite +
    draw_grob(read_panel(p$stem),
      x = tx / COMP_W, y = (ty - TITLE_H - img_h[i]) / COMP_H,
      width = col_w / COMP_W, height = img_h[i] / COMP_H, hjust = 0, vjust = 0
    ) +
    draw_label(p$tag,
      x = tx / COMP_W, y = ty / COMP_H,
      size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1
    ) +
    draw_label(p$title,
      x = (tx + 4) / COMP_W, y = ty / COMP_H,
      size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1
    ) +
    draw_label(p$subtitle,
      x = (tx + 4) / COMP_W, y = (ty - 4.2) / COMP_H,
      size = SUB_SZ, fontface = "italic", colour = "grey30", hjust = 0, vjust = 1
    )
}

ggsave(file.path(RPT_DIR, "F02_composite.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT_DIR, "F02_composite.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm", device = PDF_DEVICE, bg = "white"
)
message("F02 composite saved to ", RPT_DIR)
