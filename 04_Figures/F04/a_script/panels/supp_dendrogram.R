# F04 supplement: WGCNA gene dendrogram with the dynamic-tree-cut and merged
# module colorbands beneath it. Reads the in-memory network from setup.
pacman::p_load(here, ggplot2, png, grid)
if (!exists("wgcna_net")) source(here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))

block_genes <- wgcna_net$blockGenes[[1]]
merged_cols <- wgcna_colors[block_genes]
unmerged_cols <- wgcna_net$unmergedColors[block_genes]
color_matrix <- cbind(unmerged_cols, merged_cols)

n_mods <- length(unique(merged_cols[merged_cols != "grey"]))
n_genes <- length(merged_cols)
n_grey <- sum(merged_cols == "grey")

dendro_tmp <- tempfile(fileext = ".png")
png(dendro_tmp, width = 3200, height = 1600, res = 300)
par(mar = c(1, 4, 1, 0.5))
WGCNA::plotDendroAndColors(
  wgcna_net$dendrograms[[1]], color_matrix,
  c("Dynamic tree cut", "Merged modules"),
  main = "", dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  cex.colorLabels = 0.7, cex.axis = 0.8
)
dev.off()
dendro_img <- png::readPNG(dendro_tmp)

supp_dendro <- ggplot() +
  annotation_raster(dendro_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  labs(
    title = "WGCNA gene dendrogram and module colors",
    subtitle = sprintf(
      "Signed network; power = %d; %d modules; %s proteins (%d unassigned)",
      wgcna_power, n_mods, format(n_genes, big.mark = ","), n_grey
    )
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(face = "italic", size = 10, color = "grey30"),
    plot.margin = margin(4, 4, 4, 4)
  )

dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)
save_panel(supp_dendro, file.path(RPT_DIR, "supp", "supp_dendrogram"), 240, 120)
message("F04 supplement (dendrogram) done.")
