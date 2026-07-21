# F02 Panel A: global proteome state (PCA + PERMANOVA)
# PCA on the imputed matrix, one scatter coloured by responder group with timepoint as
# point shape and 80% group ellipses. PERMANOVA: timepoint permutes within subject
# blocks; group collapses to subject means and tests subjects as independent units, so
# the partial subjects (S28/S29) don't act as unequal whole-plot blocks. Title on the
# composite.

pacman::p_load(here, dplyr, tibble, ggplot2, vegan)

if (!exists("meta")) source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))

PA_W <- 120
PA_H <- 110

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

set.seed(42)
dist_mat <- vegdist(mat, method = "euclidean")
perm_time <- adonis2(dist_mat ~ Timepoint,
  data = pca_df, permutations = 999, strata = pca_df$Subject_ID
)
subj_mat <- rowsum(mat, pca_df$Subject_ID)
subj_mat <- subj_mat / as.integer(table(pca_df$Subject_ID)[rownames(subj_mat)])
subj_group <- pca_df$Group[match(rownames(subj_mat), pca_df$Subject_ID)]
set.seed(42)
perm_group <- adonis2(vegdist(subj_mat, method = "euclidean") ~ subj_group, permutations = 999)

sig_mark <- function(p) if (p < 0.05) "*" else "ns"
perm_label <- sprintf(
  "PERMANOVA: Group %s (p=%.2f) . Time %s (p=%.2f)",
  sig_mark(perm_group$`Pr(>F)`[1]), perm_group$`Pr(>F)`[1],
  sig_mark(perm_time$`Pr(>F)`[1]), perm_time$`Pr(>F)`[1]
)

pca_df <- pca_df |>
  mutate(gt = factor(Group_Time, levels = names(GROUP_TIME_COLORS)))

# Horizontal key above the plot: group by colour, timepoint by symbol.
xr <- range(pca_df$PC1, na.rm = TRUE)
yr <- range(pca_df$PC2, na.rm = TRUE)
ky <- yr[2] + diff(yr) * 0.14
key_grp <- tibble(
  x = xr[1] + diff(xr) * c(0.05, 0.19), lab = c("HR", "LR"),
  col = c(GROUP_COLORS[["HR"]], GROUP_COLORS[["LR"]])
)
key_tp <- tibble(
  x = xr[1] + diff(xr) * c(0.44, 0.58, 0.72), lab = c("T1", "T2", "T3"),
  shp = c(16, 17, 15)
)

pA <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color = Group, shape = Timepoint), size = 2, alpha = 0.9) +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_shape_manual(values = c(T1 = 16, T2 = 17, T3 = 15), guide = "none") +
  ggnewscale::new_scale_color() +
  stat_ellipse(aes(color = gt, group = gt), level = 0.80, linewidth = 0.45) +
  scale_color_manual(values = GROUP_TIME_COLORS, guide = "none") +
  geom_point(data = key_grp, aes(x, ky), color = key_grp$col, size = 2.6, inherit.aes = FALSE) +
  geom_text(
    data = key_grp, aes(x + diff(xr) * 0.025, ky, label = lab), hjust = 0,
    fontface = "bold", size = FIG_GEOM_TEXT - 0.4, inherit.aes = FALSE
  ) +
  geom_point(data = key_tp, aes(x, ky), shape = key_tp$shp, color = "grey30", size = 2.3, inherit.aes = FALSE) +
  geom_text(
    data = key_tp, aes(x + diff(xr) * 0.025, ky, label = lab), hjust = 0,
    size = FIG_GEOM_TEXT - 0.4, inherit.aes = FALSE
  ) +
  annotate("label",
    x = -Inf, y = -Inf, label = perm_label, hjust = -0.02, vjust = -0.05,
    size = FIG_GEOM_TEXT - 0.7, color = "grey15", fill = scales::alpha("white", 0.85),
    label.size = 0, label.padding = unit(1.5, "pt")
  ) +
  coord_cartesian(
    ylim = c(yr[1] - diff(yr) * 0.05, yr[2] + diff(yr) * 0.05), clip = "off"
  ) +
  labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
  FIG_THEME +
  theme(legend.position = "none", plot.margin = margin(t = 16, r = 5, b = 4, l = 4))

save_png(pA, file.path(RPT_DIR, "panels", "panel_a_pca"), PA_W, PA_H)
F02_AUDIT[["panel_A_pca_scores"]] <- pca_df |> select(sample, PC1, PC2, Group, Timepoint)
cat("F02 Panel A done.\n")
