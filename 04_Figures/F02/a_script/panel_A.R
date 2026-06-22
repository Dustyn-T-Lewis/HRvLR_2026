# F02 Panel A: Dual PCA (Group + Timepoint) + PERMANOVA ----
# PCA on imputed data, colored by Group then Timepoint, 80% ellipses

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
})

PA_W <- 200
PA_H <- 100
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_COLS <- GROUP_COLORS
TIME_COLS <- TIME_COLORS

# PCA on imputed matrix (samples as rows)
samp_ids <- meta$Col_ID[meta$Col_ID %in% colnames(imp_mat)]
mat <- t(imp_mat[, samp_ids])
mat <- mat[, apply(mat, 2, var, na.rm = TRUE) > 0]

pca_res <- prcomp(mat, center = TRUE, scale. = TRUE)
var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[, 1:2]) |>
  mutate(sample = rownames(pca_res$x)) |>
  left_join(meta |> select(Col_ID, Group, Timepoint, Group_Time, Subject_ID),
    by = c("sample" = "Col_ID")
  )

# PERMANOVA: Group is between-subject; Timepoint within-subject (permute within blocks)
dist_mat <- vegdist(mat, method = "euclidean")
perm_group <- adonis2(dist_mat ~ Group, data = pca_df, permutations = 999)
perm_time <- adonis2(dist_mat ~ Timepoint,
  data = pca_df, permutations = 999,
  strata = pca_df$Subject_ID
)
stat_label <- function(res, term) {
  sprintf(
    "PERMANOVA\n%s  R² = %.3f, %s", term,
    res$R2[1], fmt_p(res$`Pr(>F)`[1])
  )
}

# Shaded 80% ellipses (fill + dashed outline) for clean group separation
ellipse_layers <- function(grp) {
  list(
    stat_ellipse(aes(fill = .data[[grp]]),
      geom = "polygon",
      alpha = 0.12, level = 0.80, show.legend = FALSE
    ),
    stat_ellipse(aes(group = .data[[grp]]),
      level = 0.80,
      linewidth = 0.4, linetype = "dashed", show.legend = FALSE
    )
  )
}

# Plot 1: Group coloring
pA_group <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  ellipse_layers("Group") +
  geom_point(size = 2, alpha = 0.85) +
  annotate("text",
    x = Inf, y = Inf, label = stat_label(perm_group, "Group"),
    hjust = 1.05, vjust = 1.15, size = 2.6, color = "grey30", fontface = "bold"
  ) +
  scale_color_manual(values = GROUP_COLS) +
  scale_fill_manual(values = GROUP_COLS, guide = "none") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_pct[1]),
    y = sprintf("PC2 (%.1f%%)", var_pct[2]),
    title = "by Group", color = "Group"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

# Plot 2: Timepoint coloring
pA_time <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Timepoint)) +
  ellipse_layers("Timepoint") +
  geom_point(size = 2, alpha = 0.85) +
  annotate("text",
    x = Inf, y = Inf, label = stat_label(perm_time, "Time"),
    hjust = 1.05, vjust = 1.15, size = 2.6, color = "grey30", fontface = "bold"
  ) +
  scale_color_manual(values = TIME_COLS) +
  scale_fill_manual(values = TIME_COLS, guide = "none") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_pct[1]),
    y = sprintf("PC2 (%.1f%%)", var_pct[2]),
    title = "by Timepoint", color = "Timepoint"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

pA_combined <- (pA_group | pA_time) +
  plot_annotation(
    title = "PCA of Imputed Proteome",
    subtitle = sprintf(
      "%s proteins (imputed), %d samples | 80%% shaded ellipses",
      format(nrow(imp_df), big.mark = ","), nrow(pca_df)
    ),
    tag_levels = list(c("A", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey30"),
      plot.tag      = element_text(face = "bold", size = 15)
    )
  )

save_panel(pA_combined, file.path(RPT_DIR, "panel_a_pca"), PA_W, PA_H)
write.csv(pca_df |> select(sample, PC1, PC2, Group, Timepoint),
  file.path(DAT_DIR, "audit_panel_A_pca_scores.csv"),
  row.names = FALSE
)
cat("F02 Panel A done.\n")
