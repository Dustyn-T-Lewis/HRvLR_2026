# F02 Panel A, continuous tree: global proteome state (PCA + PERMANOVA)
# Mirrors categorical/F02_proteome/a_script/panels/panel_a_pca.R with the
# group colour/ellipse replaced by a continuous gradient on composite
# hypertrophy, and the group PERMANOVA replaced by one testing the same
# continuous phenotype (subject-collapsed, matching the group test's own
# subject-level unit) alongside timepoint. The custom hand-drawn legend key
# categorical/'s panel uses only makes sense for a small discrete key
# (HR/LR x T1/T2/T3); a continuous colour bar is a standard ggplot legend,
# so this panel keeps that instead of reproducing the manual layout.

pacman::p_load(here, dplyr, tibble, ggplot2, vegan)

if (!exists("meta")) {
  source(here(
    "03_Analysis", "continuous", "F02_proteome", "a_script", "setup.R"
  ))
}

PA_W <- 120
PA_H <- 110

samp_ids <- meta$Col_ID[meta$Col_ID %in% colnames(imp_mat)]
mat <- t(imp_mat[, samp_ids])
mat <- mat[, apply(mat, 2, var, na.rm = TRUE) > 0]

pca_res <- prcomp(mat, center = TRUE, scale. = TRUE)
var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[, 1:2]) |>
  mutate(sample = rownames(pca_res$x)) |>
  left_join(
    meta |> select(Col_ID, Timepoint, Subject_ID, comp_hypertrophy),
    by = c("sample" = "Col_ID")
  )

set.seed(42)
dist_mat <- vegdist(mat, method = "euclidean")
perm_time <- adonis2(dist_mat ~ Timepoint,
  data = pca_df, permutations = 999, strata = pca_df$Subject_ID
)
subj_mat <- rowsum(mat, pca_df$Subject_ID)
subj_mat <- subj_mat / as.integer(table(pca_df$Subject_ID)[rownames(subj_mat)])
subj_pheno <- pca_df$comp_hypertrophy[
  match(rownames(subj_mat), pca_df$Subject_ID)
]
set.seed(42)
perm_pheno <- adonis2(
  vegdist(subj_mat, method = "euclidean") ~ subj_pheno,
  permutations = 999
)

sig_mark <- function(p) if (p < 0.05) "*" else "ns"
perm_label <- sprintf(
  "PERMANOVA\nHypertrophy %s (p = %.2f)\nTime %s (p = %.2f)",
  sig_mark(perm_pheno$`Pr(>F)`[1]), perm_pheno$`Pr(>F)`[1],
  sig_mark(perm_time$`Pr(>F)`[1]), perm_time$`Pr(>F)`[1]
)

pA <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(
    aes(color = comp_hypertrophy, shape = Timepoint),
    size = 2.4, alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#B2182B", mid = "grey85", high = "#2166AC", midpoint = 0,
    name = "Composite\nhypertrophy (%)"
  ) +
  scale_shape_manual(values = c(T1 = 16, T2 = 17, T3 = 15), name = NULL) +
  annotate("label",
    x = -Inf, y = Inf, label = perm_label, hjust = -0.03, vjust = 1.08,
    size = FIG_GEOM_TEXT - 0.7, fontface = "bold", color = "grey15",
    fill = scales::alpha("white", 0.85),
    label.size = 0, label.padding = unit(2, "pt"), lineheight = 0.95
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_pct[1]),
    y = sprintf("PC2 (%.1f%%)", var_pct[2])
  ) +
  FIG_THEME +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 24, r = 5, b = 4, l = 4)
  )

save_png(pA, file.path(RPT_DIR, "panels", "panel_a_pca"), PA_W, PA_H)
F02_AUDIT[["panel_A_pca_scores"]] <- pca_df |>
  select(sample, PC1, PC2, Timepoint, comp_hypertrophy)
cat("F02 (continuous) Panel A done.\n")
