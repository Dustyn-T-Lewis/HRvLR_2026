# F02 supplement: pi-selected protein heatmaps, real arm labels beside shuffled ones
# One row per timepoint. Left column selects proteins by pi < PI_THRESH on the real arm
# contrast; right column runs the identical selection on labels shuffled across subjects.
# Rows and columns are clustered, so the blocks are the ones a reader would find. Both
# columns block, because selecting proteins for separating the arms and then displaying
# that separation is one fact shown twice. The shuffled column is the control that says so.

pacman::p_load(
  here, dplyr, tibble, limma, withr, ComplexHeatmap, circlize, grid, patchwork, ggplot2
)

if (!exists("meta")) source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))
source(here("03_Features", "contrasts.R"))
source(here("functions", "shared_pathway_utils.R"))

PI_HM_SEED <- 7
CELL_MM <- 3.4
TP_LABEL <- c(T1 = "T1 baseline", T2 = "T2 trained", T3 = "T3 acute")

goslim <- build_goslim_gene_sets(min_size = SET_FLOOR, max_size = 500)
slim_sizes <- lengths(goslim)

# A protein sits in several slim terms; take the smallest containing term so the label
# is the most specific available, and leave the rest unassigned.
slim_of <- function(genes) {
  vapply(genes, function(g) {
    hit <- names(goslim)[vapply(goslim, function(s) g %in% s, logical(1))]
    if (!length(hit)) "Unassigned" else sub("^GOSLIM_", "", hit[which.min(slim_sizes[hit])])
  }, character(1))
}

pi_selected <- function(x, g) {
  tt <- topTable(eBayes(lmFit(x, model.matrix(~g))),
    coef = 2, number = Inf, sort.by = "none"
  )
  rownames(x)[pi_score(tt$P.Value, tt$logFC) < PI_THRESH]
}

pi_panels <- lapply(names(TP_LABEL), function(tp) {
  samples <- meta$Col_ID[meta$Timepoint == tp]
  x <- imp_mat[, samples]
  real <- as.character(meta$Group[match(samples, meta$Col_ID)])
  shuffled <- with_seed(PI_HM_SEED + match(tp, names(TP_LABEL)), sample(real))

  lapply(
    list(
      list(tag = "Selected on arm", lab = real),
      list(tag = "Selected on shuffle", lab = shuffled)
    ),
    function(s) {
      keep <- pi_selected(x, factor(s$lab))
      z <- t(scale(t(x[keep, , drop = FALSE])))
      colnames(z) <- samples
      list(
        z = z, arm = s$lab, slim = slim_of(rownames(z)),
        title = sprintf("%s  |  %s  |  %d proteins", TP_LABEL[[tp]], s$tag, nrow(z))
      )
    }
  )
})
pi_panels <- unlist(pi_panels, recursive = FALSE)

# One colour per GO Slim term seen anywhere, so the same term reads the same in
# every panel.
slim_levels <- sort(unique(unlist(lapply(pi_panels, function(p) p$slim))))
slim_cols <- stats::setNames(
  grDevices::hcl.colors(length(slim_levels), "Dark 3"), slim_levels
)

z_scale <- colorRamp2(
  c(-2, 0, 2), c(DIR_COLORS[["Down"]], "white", DIR_COLORS[["Up"]])
)

draw_pi_panel <- function(p) {
  ht <- Heatmap(
    p$z,
    name = "z", col = z_scale,
    column_title = p$title, column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_dend = TRUE, show_column_dend = TRUE,
    row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"),
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 5.6), column_names_gp = gpar(fontsize = 5.6),
    width = unit(ncol(p$z) * CELL_MM, "mm"),
    height = unit(nrow(p$z) * CELL_MM, "mm"),
    top_annotation = HeatmapAnnotation(
      Arm = p$arm,
      col = list(Arm = GROUP_COLORS),
      annotation_name_gp = gpar(fontsize = 6),
      simple_anno_size = unit(2.6, "mm")
    ),
    left_annotation = rowAnnotation(
      `GO Slim` = p$slim,
      col = list(`GO Slim` = slim_cols),
      annotation_name_gp = gpar(fontsize = 6),
      simple_anno_size = unit(2.6, "mm"),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.2),
        grid_width = unit(2.5, "mm"), grid_height = unit(2.5, "mm")
      )
    ),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.6),
      grid_width = unit(2.5, "mm"), legend_height = unit(14, "mm")
    )
  )
  wrap_elements(grid.grabExpr(draw(ht, merge_legend = TRUE)))
}

# Rows sized by the taller panel in each pair so cells stay square across the sheet.
row_heights <- vapply(seq(1, length(pi_panels), by = 2), function(i) {
  max(nrow(pi_panels[[i]]$z), nrow(pi_panels[[i + 1]]$z))
}, numeric(1))

p_pi <- wrap_plots(lapply(pi_panels, draw_pi_panel), ncol = 2) +
  plot_layout(heights = row_heights) +
  plot_annotation(
    title = "pi-selected proteins, real arm labels beside shuffled ones",
    caption = paste(
      "Left selects proteins by pi < 0.05 on the real arm contrast; right runs the",
      "identical selection on labels shuffled across subjects.\nRows and columns are",
      "clustered within each panel. Both block, and at every timepoint the shuffle",
      "selects more proteins than the arm does."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 9),
      plot.caption = element_text(hjust = 0, size = 6, colour = "grey35")
    )
  )

# Each panel needs its cells plus about 21 mm of dendrogram, title and column labels.
save_png(
  p_pi, file.path(RPT_DIR, "supp", "supp_pi_heatmap"),
  200, sum(row_heights) * CELL_MM + 21 * length(row_heights)
)
F02_AUDIT[["supp_pi_heatmap"]] <- bind_rows(lapply(pi_panels, function(p) {
  tibble(panel = p$title, gene = rownames(p$z), slim = p$slim)
}))
cat("F02 pi heatmap supplement done.\n")
