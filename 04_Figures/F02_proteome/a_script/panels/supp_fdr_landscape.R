# F02 supplement: every protein, every contrast, nothing selected
# All 1900 proteins by the nine contrasts, coloured by logFC and clustered on rows.
# No feature is chosen on the outcome, so there is no selection to control for and no
# shuffled counterpart is needed. The strip on the right is each protein's smallest BH
# q across the nine; it never reaches the threshold, which is the result the sheet exists
# to show.

pacman::p_load(
  here, dplyr, readr, ComplexHeatmap, circlize, grid, patchwork, ggplot2
)

if (!exists("meta")) source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))

FDR_CONTRASTS <- c(
  Baseline_HRvLR = "HR - LR", Trained_HRvLR = "HR - LR", Acute_HRvLR = "HR - LR",
  Training_Interaction = "Interaction", Acute_Interaction = "Interaction",
  Training_HR = "Within arm", Training_LR = "Within arm",
  Acute_HR = "Within arm", Acute_LR = "Within arm"
)

lfc <- as.matrix(dep_df[, paste0("logFC_", names(FDR_CONTRASTS))])
qval <- as.matrix(dep_df[, paste0("adj.P.Val_", names(FDR_CONTRASTS))])
colnames(lfc) <- names(FDR_CONTRASTS)
rownames(lfc) <- dep_df$gene

keep <- stats::complete.cases(lfc)
lfc <- lfc[keep, ]
min_q <- apply(qval[keep, ], 1, min, na.rm = TRUE)

lim <- stats::quantile(abs(lfc), 0.99)
ht <- Heatmap(
  lfc,
  name = "logFC",
  col = colorRamp2(c(-lim, 0, lim), c(DIR_COLORS[["Down"]], "white", DIR_COLORS[["Up"]])),
  column_split = factor(
    unname(FDR_CONTRASTS),
    levels = c("HR - LR", "Interaction", "Within arm")
  ),
  cluster_rows = TRUE, cluster_columns = FALSE,
  show_row_names = FALSE, show_row_dend = FALSE,
  column_names_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 7, fontface = "bold"),
  # cairo is absent on this machine, so the default raster device cannot open its
  # temp png; ragg is present and writes the same thing.
  use_raster = TRUE, raster_quality = 3, raster_device = "agg_png",
  width = unit(length(FDR_CONTRASTS) * 7, "mm"), height = unit(105, "mm"),
  right_annotation = rowAnnotation(
    `min q` = min_q,
    col = list(`min q` = colorRamp2(c(0, 0.05, 1), c("#B2182B", "#F4A582", "grey92"))),
    annotation_name_gp = gpar(fontsize = 6),
    simple_anno_size = unit(3.5, "mm"),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.6),
      grid_width = unit(2.5, "mm"), legend_height = unit(20, "mm")
    )
  ),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.6),
    grid_width = unit(2.5, "mm"), legend_height = unit(20, "mm")
  )
)

p_fdr <- wrap_elements(grid.grabExpr(draw(ht, merge_legend = TRUE))) +
  plot_annotation(
    caption = sprintf(
      paste(
        "All %s proteins, no selection, rows clustered on the nine logFC profiles.",
        "Smallest BH q anywhere is %.3f, in %s;\nnothing crosses 0.05, so the q strip",
        "stays grey. Colour is capped at the 99th percentile of |logFC|."
      ),
      format(nrow(lfc), big.mark = ","), min(min_q),
      names(FDR_CONTRASTS)[which.min(apply(qval[keep, ], 2, min, na.rm = TRUE))]
    ),
    theme = theme(plot.caption = element_text(hjust = 0, size = 6, colour = "grey35"))
  )

save_png(p_fdr, file.path(RPT_DIR, "supp", "supp_fdr_landscape"), 150, 135)
F02_AUDIT[["supp_fdr_landscape"]] <- tibble::tibble(
  gene = rownames(lfc), min_q = min_q
) |>
  arrange(min_q)
cat("F02 FDR landscape supplement done.\n")
