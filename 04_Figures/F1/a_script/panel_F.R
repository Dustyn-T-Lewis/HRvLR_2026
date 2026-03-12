# Figure 1 — Panel F: fGSEA Grouped Bar Chart (Pathway Enrichment) ----

if (!exists("meta")) source("04_Figures/F1/a_script/HRvLR_F1_setup.R")
source("04_Figures/shared/pathway_utils.R")

RPT <- RPT_DIR
DAT <- DAT_DIR

PF_W <- 120

# Build pathway collection ----

pw_collection <- build_pathway_collection(min_size = 10, max_size = 500)

# Run fGSEA per contrast ----

set.seed(42)
fgsea_all <- list()

for (ctr in CONTRASTS) {
  stats <- dep_df[[paste0("t_", ctr)]]
  names(stats) <- dep_df$gene
  stats <- sort(stats[!is.na(stats) & is.finite(stats)], decreasing = TRUE)

  res <- run_fgsea_deduplicated(
    ranks          = stats,
    pathways       = pw_collection,
    jaccard_cutoff = 0.5,
    nperm          = 10000,
    min_size       = 10,
    max_size       = 500
  )
  res$contrast <- ctr
  fgsea_all[[ctr]] <- res
  message(sprintf("  fgsea done: %s", ctr))
}

fgsea_combined <- dplyr::bind_rows(fgsea_all)

# Export full results ----

fgsea_export <- fgsea_combined |>
  mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) |>
  arrange(database, contrast, padj)
write_csv(fgsea_export, file.path(DAT, "06_panel_F_fgsea_results.csv"))

# Count significant up/down per database and contrast ----

DISPLAY_DBS <- c("Hallmark", "KEGG", "Reactome", "WikiPathways", "BioCarta", "PID", "GO:BP")

count_df <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05, database %in% DISPLAY_DBS) |>
  group_by(contrast, database) |>
  summarise(
    Up   = sum(NES > 0),
    Down = sum(NES < 0),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(cols = c(Up, Down), names_to = "direction",
                      values_to = "count")

nonempty_dbs <- count_df |>
  group_by(database) |> filter(sum(count) > 0) |> pull(database) |> unique()
count_df <- count_df |> filter(database %in% nonempty_dbs)

count_df$contrast  <- factor(count_df$contrast, levels = DISPLAY_ORDER)
count_df$database  <- factor(count_df$database, levels = intersect(DISPLAY_DBS, nonempty_dbs))
count_df$direction <- factor(count_df$direction, levels = c("Up", "Down"))

n_facets <- length(levels(count_df$database))
PF_H <- max(100, n_facets * 40)

lbl_sz <- scale_text(BASE_COUNT, PF_W)

# Background rects ----
n_ctr <- length(DISPLAY_ORDER)
bg_rects <- lapply(seq_len(n_ctr), function(i) {
  annotate("rect",
           xmin = i - 0.5, xmax = i + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS[DISPLAY_ORDER[i]], alpha = 0.20,
           color = "grey70", linewidth = 0.2)
})

# Plot ----

pF <- ggplot(count_df, aes(x = contrast, y = count, fill = direction))
for (r in bg_rects) pF <- pF + r
pF <- pF +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(data = \(d) d |> filter(count >= 5) |> mutate(mid = count / 2),
            aes(y = mid, label = count),
            position = position_dodge(width = 0.7),
            vjust = 0.5, hjust = 0.5, size = lbl_sz,
            color = "white", fontface = "bold") +
  geom_text(data = \(d) d |> filter(count > 0, count < 5),
            aes(y = count, label = count),
            position = position_dodge(width = 0.7),
            vjust = -0.3, hjust = 0.5, size = lbl_sz,
            color = "grey30", fontface = "bold") +
  facet_grid(database ~ ., scales = "free_y") +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "Pathway Enrichment",
       subtitle = "fGSEA (padj < 0.05); Jaccard\ndedup (threshold 0.5)",
       x = NULL, y = "Significant pathways",
       tag = "F") +
  FIG_THEME +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 6.5),
        legend.position = "none",
        strip.text.y   = element_text(size = 7, face = "bold", angle = 0))

# Export NES summary ----

nes_summary <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05) |>
  group_by(contrast, database) |>
  summarise(
    n_sig = n(), n_up = sum(NES > 0), n_down = sum(NES < 0),
    median_NES = median(NES), mean_NES = mean(NES), sd_NES = sd(NES),
    min_padj = min(padj), median_padj = median(padj),
    .groups = "drop"
  )
write_csv(nes_summary, file.path(DAT, "audit_panel_F_nes_summary.csv"))

save_panel(pF, file.path(RPT, "panel_F_fgsea"),
           width = PF_W, height = PF_H)
