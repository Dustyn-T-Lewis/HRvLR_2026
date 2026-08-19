# F02 Panel H: unsupervised sample clustering on the full proteome
# Sample-sample correlation over all 1900 proteins with no feature
# selection, so the arms cannot be separated by construction. Rows and
# columns are ordered by hierarchical clustering; the strips above carry
# arm and timepoint. If responder status organised the proteome the blue
# and red strip would follow the dendrogram, and it does not.
# Title on the composite.

pacman::p_load(here, dplyr, tibble, tidyr, ggplot2, patchwork)

if (!exists("meta")) {
  source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))
}

PH_W <- 130
PH_H <- 110

cor_mat <- cor(imp_mat, method = "pearson")
ord <- hclust(as.dist(1 - cor_mat), method = "average")$order
samp_order <- colnames(cor_mat)[ord]

ann <- meta |>
  filter(Col_ID %in% samp_order) |>
  mutate(Col_ID = factor(Col_ID, levels = samp_order))

# Does a sample's nearest neighbour share its arm? Same-subject pairs are
# masked first: a subject's other timepoints are its closest samples by
# a wide margin and they always share its arm, so leaving them in
# measures subject identity, not arm.
arm <- stats::setNames(as.character(meta$Group), meta$Col_ID)
subj <- stats::setNames(as.character(meta$subject), meta$Col_ID)

nn_mat <- cor_mat
nn_mat[outer(subj[rownames(nn_mat)], subj[colnames(nn_mat)], "==")] <- -Inf
nn <- rownames(nn_mat)[apply(nn_mat, 2, which.max)]
nn_share <- mean(arm[nn] == arm[colnames(nn_mat)])

# Chance is the share of eligible (different-subject) partners in the same arm,
# averaged over samples, not a naive marginal.
expected <- mean(vapply(colnames(nn_mat), function(s) {
  elig <- subj != subj[[s]]
  mean(arm[elig] == arm[[s]])
}, numeric(1)))

heat <- as.data.frame(as.table(cor(imp_mat))) |>
  rlang::set_names(c("row", "col", "r")) |>
  mutate(
    row = factor(row, levels = samp_order),
    col = factor(col, levels = samp_order)
  )

p_heat <- ggplot(heat, aes(col, row, fill = r)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c("#4393C3", "white", "#B2182B"),
    limits = range(heat$r), name = "Pearson r"
  ) +
  labs(
    x = NULL, y = NULL,
    caption = sprintf(
      paste(
        "Nearest neighbour in another subject shares arm for %d of %d samples,",
        "where chance predicts %.0f."
      ),
      round(nn_share * ncol(nn_mat)), ncol(nn_mat), expected * ncol(nn_mat)
    )
  ) +
  FIG_THEME +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), legend.position = "right",
    plot.caption = element_text(
      hjust = 0, size = FIG_GEOM_TEXT - 0.4, colour = "grey35"
    )
  )

strip <- function(fill_var, values, name) {
  ggplot(ann, aes(Col_ID, 1, fill = .data[[fill_var]])) +
    geom_raster() +
    scale_fill_manual(values = values, name = name) +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), legend.position = "right",
      plot.margin = margin(0, 0, 1, 0)
    )
}

p_h <- strip("Group", GROUP_COLORS, "Arm") /
  strip("Timepoint", TIME_COLORS, "Timepoint") /
  p_heat +
  plot_layout(heights = c(1, 1, 26), guides = "collect")

save_png(p_h, file.path(RPT_DIR, "panels", "panel_h_clustering"), PH_W, PH_H)
F02_AUDIT[["panel_H_sample_correlation"]] <- tibble(
  sample = colnames(cor_mat), nearest_neighbour = nn,
  arm = arm[colnames(cor_mat)], nn_arm = arm[nn]
)
cat("F02 Panel H done.\n")
