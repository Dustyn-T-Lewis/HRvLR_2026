# Figure 1 — Panel C: Dual PCA (Group + Timepoint) + PERMANOVA ----

if (!exists("meta")) source("04_Figures/F1/a_script/HRvLR_F1_setup.R")

PC_W <- 200
PC_H <- 125

# ── Computation ----

pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

set.seed(42)
n_boot <- 1000
n_prot <- nrow(imp_mat)
boot_var <- matrix(NA_real_, nrow = n_boot, ncol = 2)
for (b in seq_len(n_boot)) {
  idx <- sample(n_prot, replace = TRUE)
  bp  <- prcomp(t(imp_mat[idx, ]), center = TRUE, scale. = TRUE)
  boot_var[b, ] <- 100 * summary(bp)$importance[2, 1:2]
}
var_ci <- data.frame(
  PC      = c("PC1", "PC2"),
  var_pct = var_pct,
  ci_lo   = apply(boot_var, 2, quantile, 0.025),
  ci_hi   = apply(boot_var, 2, quantile, 0.975)
)

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(sample_id = rownames(pca$x)) |>
  left_join(meta, by = "sample_id")

dist_mat <- dist(scale(t(imp_mat)))
set.seed(42)
perm_res <- adonis2(dist_mat ~ Group * Timepoint, data = meta,
                    permutations = how(nperm = 999, blocks = meta$subject),
                    by = "terms")

perm_terms <- c("Group", "Timepoint", "Group:Timepoint")
perm_r2 <- perm_res[perm_terms, "R2"]
perm_pv <- perm_res[perm_terms, "Pr(>F)"]
.bare_p <- function(p) {
  if (is.na(p)) return("NA")
  sub("^p [=<] ", "", fmt_p(p))
}

bd_grp  <- betadisper(dist_mat, meta$Group)
bd_time <- betadisper(dist_mat, meta$Timepoint)
bd_gt   <- betadisper(dist_mat, meta$group)
bd_grp_p  <- permutest(bd_grp,  pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
bd_time_p <- permutest(bd_time, pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
bd_gt_p   <- permutest(bd_gt,   pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
if (bd_grp_p < 0.05 || bd_time_p < 0.05)
  warning("Heterogeneous dispersions detected — interpret PERMANOVA with caution")

# Shared axis labels
x_lab <- sprintf("PC1 (%.1f%% [%.1f, %.1f])", var_pct[1], var_ci$ci_lo[1], var_ci$ci_hi[1])
y_lab <- sprintf("PC2 (%.1f%% [%.1f, %.1f])", var_pct[2], var_ci$ci_lo[2], var_ci$ci_hi[2])

# ── Left PCA: colored by Group ----

perm_label_grp <- sprintf(
  "PERMANOVA\nGroup  R\u00b2=%.3f, p=%s\nG\u00d7T  R\u00b2=%.3f, p=%s",
  perm_r2[1], .bare_p(perm_pv[1]),
  perm_r2[3], .bare_p(perm_pv[3]))

pC_group <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  stat_ellipse(aes(fill = Group), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = Group), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(aes(shape = Timepoint), size = 2.0, alpha = 0.85) +
  annotate("text", x = -Inf, y = Inf, label = perm_label_grp,
           hjust = -0.05, vjust = 1.15,
           size = scale_text(BASE_STAT - 1.2, PC_W / 2), color = "grey30",
           fontface = "bold") +
  scale_color_manual(values = GROUP_COLORS) +
  scale_fill_manual(values = GROUP_COLORS) +
  scale_shape_manual(values = SHAPE_TP) +
  labs(x = x_lab, y = y_lab) +
  FIG_THEME +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6.5))

# ── Right PCA: colored by Timepoint ----

TP_COLORS <- c(T1 = "#999999", T2 = "#E69F00", T3 = "#8B5CF6")

perm_label_time <- sprintf(
  "PERMANOVA\nTime  R\u00b2=%.3f, p=%s\nG\u00d7T  R\u00b2=%.3f, p=%s",
  perm_r2[2], .bare_p(perm_pv[2]),
  perm_r2[3], .bare_p(perm_pv[3]))

pC_time <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Timepoint)) +
  stat_ellipse(aes(fill = Timepoint), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = Timepoint), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(aes(shape = Group), size = 2.0, alpha = 0.85) +
  annotate("text", x = -Inf, y = Inf, label = perm_label_time,
           hjust = -0.05, vjust = 1.15,
           size = scale_text(BASE_STAT - 1.2, PC_W / 2), color = "grey30",
           fontface = "bold") +
  scale_color_manual(values = TP_COLORS) +
  scale_fill_manual(values = TP_COLORS) +
  scale_shape_manual(values = c(HR = 16, LR = 17)) +
  labs(x = x_lab, y = y_lab) +
  FIG_THEME +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6.5))

# ── Compose side-by-side ----

pC_combined <- (pC_group | pC_time) +
  plot_annotation(
    title = "Principal Component Analysis (PCA)",
    subtitle = sprintf("%s proteins (HPA-filtered, cycloess, %s-imputed); n = %d samples",
                       format(nrow(imp_df), big.mark = ","), BEST_IMP_METHOD, nrow(meta)),
    tag_levels = list(c("C", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                   color = "grey30"),
      plot.tag      = element_text(face = "bold", size = FIG_TAG_SIZE)
    )
  )

# ── Export ----

write_csv(var_ci, file.path(DAT_DIR, "audit_panel_C_pca_variance_ci.csv"))
betadisper_results <- data.frame(
  factor      = c("Group", "Timepoint", "Group_Time"),
  p_value     = c(bd_grp_p, bd_time_p, bd_gt_p),
  significant = c(bd_grp_p < 0.05, bd_time_p < 0.05, bd_gt_p < 0.05)
)
write_csv(betadisper_results, file.path(DAT_DIR, "audit_panel_C_betadisper.csv"))

save_panel(pC_combined, file.path(RPT_DIR, "panel_C_pca"),
           width = PC_W, height = PC_H)
