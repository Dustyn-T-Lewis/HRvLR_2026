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

PA_W <- 200; PA_H <- 100
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_COLS <- c(HR = "#2166AC", LR = "#B2182B")
TIME_COLS  <- c(T1 = "#F4A582", T2 = "#4393C3")

# PCA on imputed matrix (samples as rows)
samp_ids <- meta$Col_ID[meta$Col_ID %in% colnames(imp_mat)]
mat <- t(imp_mat[, samp_ids])
mat <- mat[, apply(mat, 2, var, na.rm = TRUE) > 0]

pca_res <- prcomp(mat, center = TRUE, scale. = TRUE)
var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[, 1:2]) |>
  mutate(sample = rownames(pca_res$x)) |>
  left_join(meta |> select(Col_ID, Group, Timepoint, Group_Time, Subject_ID),
            by = c("sample" = "Col_ID"))

# PERMANOVA
dist_mat <- vegdist(mat, method = "euclidean")
perm_res <- adonis2(dist_mat ~ Group * Timepoint,
                    data = pca_df, permutations = 999,
                    strata = pca_df$Subject_ID)
perm_txt <- sprintf("PERMANOVA: Group R²=%.3f %s | Time R²=%.3f %s",
                    perm_res$R2[1], ifelse(perm_res$`Pr(>F)`[1] < 0.05, "*", "ns"),
                    perm_res$R2[2], ifelse(perm_res$`Pr(>F)`[2] < 0.05, "*", "ns"))

# Plot 1: Group coloring
pA_group <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  stat_ellipse(level = 0.80, linewidth = 0.5, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = GROUP_COLS) +
  labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", var_pct[2]),
       title = "by Group", color = "Group") +
  FIG_THEME +
  theme(legend.position = "bottom")

# Plot 2: Timepoint coloring
pA_time <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Timepoint)) +
  stat_ellipse(level = 0.80, linewidth = 0.5, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = TIME_COLS) +
  labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]),
       y = sprintf("PC2 (%.1f%%)", var_pct[2]),
       title = "by Timepoint", color = "Timepoint") +
  FIG_THEME +
  theme(legend.position = "bottom")

pA_combined <- (pA_group | pA_time) +
  plot_annotation(
    title    = "PCA of Imputed Proteome",
    subtitle = perm_txt,
    tag_levels = list(c("A", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "grey30"),
      plot.tag      = element_text(face = "bold", size = 15)
    )
  )

save_panel(pA_combined, file.path(RPT_DIR, "panel_A_pca_MAIN"), PA_W, PA_H)
write.csv(pca_df |> select(sample, PC1, PC2, Group, Timepoint),
          file.path(DAT_DIR, "audit_panel_A_pca_scores.csv"), row.names = FALSE)
cat("F02 Panel A done.\n")
