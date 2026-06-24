# F03 Panel C: Mfuzz gap-trajectory structure ----
setwd(rprojroot::find_rstudio_root_file())
if (!exists("mfuzz_mem")) source("04_Figures/F03/a_script/HRvLR_F03_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# Normalise group_id (bare integers 1..4) to c1..c4 for CLUSTER_COLORS
mfuzz_mem <- mfuzz_mem |>
  mutate(cluster = paste0("c", group_id))

# Z-score each protein across its 3 timepoints (the space Mfuzz clustered on)
sub_mat <- gap_mat[mfuzz_mem$protein_id, , drop = FALSE]
gap_z <- t(scale(t(sub_mat)))

# PCA on z-scored gap trajectories
pca <- prcomp(gap_z, center = FALSE, scale. = FALSE)
var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(protein_id = rownames(pca$x)) |>
  left_join(
    mfuzz_mem |> select(protein_id, cluster, membership = membership_weight),
    by = "protein_id"
  )

# Cluster centroids in PCA space
pca_centroids <- pca_df |>
  group_by(cluster) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

# Cluster centroid trajectories over T1/T2/T3
centroids <- as.data.frame(gap_z) |>
  mutate(protein_id = rownames(gap_z)) |>
  left_join(mfuzz_mem |> select(protein_id, cluster), by = "protein_id") |>
  group_by(cluster) |>
  summarise(T1 = mean(T1), T2 = mean(T2), T3 = mean(T3), .groups = "drop") |>
  pivot_longer(c(T1, T2, T3), names_to = "timepoint", values_to = "mean_z") |>
  mutate(timepoint = factor(timepoint, levels = c("T1", "T2", "T3")))

# Left: PCA shaded by fuzzy cluster
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = cluster, alpha = membership)) +
  geom_point(size = 1.6) +
  geom_text(
    data = pca_centroids,
    aes(PC1, PC2, label = cluster, colour = cluster),
    inherit.aes = FALSE,
    fontface = "bold", size = FIG_GEOM_TEXT + 0.8
  ) +
  scale_colour_manual(values = CLUSTER_COLORS, name = "Cluster") +
  scale_alpha(range = c(0.25, 1), name = "Membership") +
  labs(
    x = sprintf("PC1 (%.0f%%)", var_pct[1]),
    y = sprintf("PC2 (%.0f%%)", var_pct[2])
  ) +
  FIG_THEME

# Right: cluster centroid trajectories
p_traj <- ggplot(centroids, aes(timepoint, mean_z, colour = cluster, group = cluster)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_colour_manual(values = CLUSTER_COLORS, guide = "none") +
  labs(x = NULL, y = "HR-LR gap (z, per protein)") +
  FIG_THEME

panel_c <- (p_pca | p_traj) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "Temporal clustering of the HR-LR proteome gap",
    subtitle = sprintf(
      "Fuzzy c-means on the T1/T2/T3 gap trajectory; point opacity = cluster membership (n approx 16; descriptive)"
    ),
    tag_levels = list(c("C", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE, color = "grey30"),
      plot.tag      = element_text(face = "bold", size = FIG_TAG_SIZE)
    )
  )

save_panel(panel_c, file.path(RPT_DIR, "panel_c_mfuzz_pca"), 200, 95)
write.csv(centroids, file.path(DAT_DIR, "audit_panel_C_cluster_centroids.csv"), row.names = FALSE)
cat("F03 Panel C done.\n")
