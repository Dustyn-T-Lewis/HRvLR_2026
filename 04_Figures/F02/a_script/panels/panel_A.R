# F02 Panel A: global proteome state (PCA + PERMANOVA)
# PCA on the imputed matrix, one scatter coloured by responder group with timepoint as
# point shape and 80% group ellipses. PERMANOVA: timepoint permutes within subject
# blocks; group collapses to subject means and tests subjects as independent units, so
# the partial subjects (S28/S29) don't act as unequal whole-plot blocks. Title on the
# composite.

pacman::p_load(here, dplyr, tibble, ggplot2, vegan)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

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
perm_group <- adonis2(vegdist(subj_mat, method = "euclidean") ~ subj_group, permutations = 999)

perm_row <- function(res, term) {
  sprintf("%-6s R²=%.3f  F=%.2f  %s", term, res$R2[1], res$F[1], fmt_p(res$`Pr(>F)`[1]))
}
perm_label <- paste(
  "PERMANOVA (999 perm.)",
  perm_row(perm_group, "Group"),
  perm_row(perm_time, "Time"),
  sep = "\n"
)

# In-panel symbol key placed in a blank right strip so it never overlaps points.
xr <- range(pca_df$PC1, na.rm = TRUE)
yr <- range(pca_df$PC2, na.rm = TRUE)
kx <- xr[2] + 0.07 * diff(xr)
lx <- kx + 0.05 * diff(xr)
yk <- yr[2] - diff(yr) * c(0.00, 0.11, 0.21, 0.37, 0.48, 0.58, 0.68)
key_grp <- tibble(x = kx, lx = lx, y = yk[c(2, 3)], lab = c("HR", "LR"))
key_tp <- tibble(x = kx, lx = lx, y = yk[c(5, 6, 7)], lab = c("T1", "T2", "T3"))
key_hdr <- tibble(x = kx - 0.015 * diff(xr), y = yk[c(1, 4)], lab = c("Group", "Time"))

pA <- ggplot(pca_df, aes(PC1, PC2)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.10, level = 0.80, color = NA) +
  stat_ellipse(aes(color = Group, group = Group), level = 0.80, linewidth = 0.4, linetype = "dashed") +
  geom_point(aes(color = Group, shape = Timepoint), size = 2, alpha = 0.85) +
  geom_point(data = key_grp, aes(x, y, color = lab), shape = 16, size = 2, inherit.aes = FALSE) +
  geom_point(data = key_tp, aes(x, y, shape = lab), color = "grey35", size = 2, inherit.aes = FALSE) +
  geom_text(
    data = rbind(key_grp[c("lx", "y", "lab")], key_tp[c("lx", "y", "lab")]),
    aes(lx, y, label = lab), hjust = 0, size = FIG_GEOM_TEXT - 0.5, inherit.aes = FALSE
  ) +
  geom_text(data = key_hdr, aes(x, y, label = lab), hjust = 0, fontface = "bold", size = FIG_GEOM_TEXT - 0.3, inherit.aes = FALSE) +
  annotate("label",
    x = -Inf, y = Inf, label = perm_label, hjust = -0.015, vjust = 1.02,
    size = FIG_GEOM_TEXT - 1, color = "black", family = "mono", lineheight = 1.05,
    fill = scales::alpha("white", 0.85), label.size = 0, label.padding = unit(1, "pt")
  ) +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_shape_manual(values = c(T1 = 16, T2 = 17, T3 = 15), guide = "none") +
  coord_cartesian(xlim = c(xr[1], lx + 0.12 * diff(xr)), clip = "off") +
  labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
  FIG_THEME +
  theme(legend.position = "none")

save_png(pA, file.path(RPT_DIR, "panels", "panel_a_pca"), PA_W, PA_H)
write.csv(pca_df |> select(sample, PC1, PC2, Group, Timepoint),
  file.path(DAT_DIR, "audit_panel_A_pca_scores.csv"),
  row.names = FALSE
)
cat("F02 Panel A done.\n")
