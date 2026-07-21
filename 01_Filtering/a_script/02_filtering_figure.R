#!/usr/bin/env Rscript
# Filtering QC figure. Reports what was removed and what contamination remained, rather than
# deleting the evidence and asserting the samples were clean.

pacman::p_load(here, readr, dplyr, tidyr, stringr, ggplot2, forcats, patchwork, scales)

# geom_jitter draws from the RNG, so without a seed this figure was a different PNG every run
# and could never be checked for drift. The seed fixes only where the jittered points land.
set.seed(42)

fx <- readRDS(here("01_Filtering", "c_data", "filtering_intermediates.rds"))
report_dir <- here("01_Filtering", "b_reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

waterfall <- fx$filter_log |>
  filter(!is.na(n_removed)) |>
  mutate(
    step = fct_inorder(str_remove(step, "^remove: ")),
    ymin = n_after,
    ymax = n_after + n_removed
  )

panel_a <- ggplot(waterfall, aes(step)) +
  geom_rect(aes(xmin = as.numeric(step) - 0.35, xmax = as.numeric(step) + 0.35, ymin = ymin, ymax = ymax),
    fill = "#D6604D", alpha = 0.9
  ) +
  geom_step(aes(y = n_after, group = 1), direction = "mid", colour = "grey35", linewidth = 0.5) +
  geom_point(aes(y = n_after), colour = "grey25", size = 1.6) +
  geom_text(aes(y = ymax, label = if_else(n_removed > 0, paste0("-", n_removed), "")),
    vjust = -0.6, size = 3, colour = "#B2182B"
  ) +
  scale_y_continuous(labels = comma, limits = c(1800, NA)) +
  labs(
    title = "A  Filtering cascade",
    subtitle = sprintf(
      "%s proteins measured -> %s analysed. Red = removed at that step.",
      comma(fx$n_raw), comma(tail(fx$filter_log$n_after, 1))
    ),
    x = NULL, y = "proteins retained"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

qc <- fx$qc_index |>
  filter(panel != "myofibre") |>
  mutate(
    panel = factor(panel, levels = c("blood", "leukocyte", "keratin", "adipose")),
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  )

panel_b <- ggplot(qc, aes(Timepoint, pct_signal, fill = Timepoint)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
  geom_jitter(aes(shape = Group), width = 0.15, size = 1.4, alpha = 0.85) +
  facet_wrap(~panel, scales = "free_y", nrow = 1) +
  scale_fill_brewer(palette = "Blues", guide = "none") +
  scale_shape_manual(values = c(HR = 16, LR = 1), name = NULL) +
  labs(
    title = "B  Contamination per sample, measured before removal",
    subtitle = paste(
      "Blood rises at T3 (hyperaemia): +15.3 pp from T2, 14/16 subjects, paired p = 0.0006, and leukocytes rise with it",
      "(post-exercise infiltration). Both rises are shared by the arms (group x timepoint p = 0.25), so they cancel in the",
      "interaction contrast but confound the per-arm acute contrasts. Keratin and adipose are negligible throughout.",
      sep = "\n"
    ),
    x = NULL, y = "% of sample signal"
  )

flags <- fx$outlier_diag |>
  filter(n_flags > 0) |>
  select(Col_ID, miss_flag, pca_flag, mad_flag, cor_flag, n_flags, consensus_outlier) |>
  pivot_longer(ends_with("_flag"), names_to = "method", values_to = "flagged") |>
  mutate(
    method = recode(method,
      miss_flag = "missingness", pca_flag = "PCA distance",
      mad_flag = "median MAD", cor_flag = "inter-sample cor"
    ),
    Col_ID = fct_reorder(Col_ID, n_flags)
  )

panel_c <- ggplot(flags, aes(method, Col_ID, fill = flagged)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    data = distinct(flags, Col_ID, consensus_outlier),
    aes(x = 4.6, y = Col_ID, label = if_else(consensus_outlier, "dropped", "kept"), colour = consensus_outlier),
    inherit.aes = FALSE, hjust = 0, size = 3, show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(`TRUE` = "#D6604D", `FALSE` = "grey92"), name = NULL,
    labels = c(`TRUE` = "flagged", `FALSE` = "clean")
  ) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey45")) +
  coord_cartesian(xlim = c(0.5, 6.4), clip = "off") +
  labs(
    title = "C  Outlier consensus, flagged samples only",
    subtitle = sprintf(
      "%d of 48 samples drew any flag; %d met the >=%d-of-4 consensus and were dropped. All three are HR.",
      n_distinct(flags$Col_ID), length(fx$outlier_ids), fx$cfg$outlier_k
    ),
    x = NULL, y = NULL
  )

fig <- panel_a / panel_b / panel_c +
  plot_layout(heights = c(1, 1.2, 0.9)) &
  theme_minimal(base_size = 10) &
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey35"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(report_dir, "filtering_qc.pdf"), fig, width = 10, height = 11)
ggsave(file.path(report_dir, "filtering_qc.png"), fig, width = 10, height = 11, dpi = 200)
cat("wrote", file.path(report_dir, "filtering_qc.png"), "\n")

# Supplement: why the data-driven blood-index cut beats the erythrocyte-RNA cut. Each protein's
# correlation with the per-sample blood index (x) against HPA's erythrocyte/myonuclei RNA ratio
# (y). The blood cut cleanly separates red-cell from muscle where the RNA ratio traps ribosomes.
muscle_markers <- c(
  "MB", "CKM", "ACTA1", "MYH1", "MYH2", "MYH7", "ALDOA", "CASQ1", "PYGM",
  "ATP2A1", "DES", "TNNT3", "TNNI2", "MYBPC1", "TNNC2", "MYL1"
)
bi_df <- fx$protein_calls |>
  filter(!is.na(blood_cor)) |>
  mutate(
    class = case_when(
      verdict == "remove: blood-tracking" ~ "removed: blood-tracking",
      gene %in% muscle_markers ~ "muscle marker",
      str_detect(gene, "^RP[LS]") ~ "ribosomal",
      TRUE ~ "other"
    ),
    ratio = log2((ery + 1) / (myo + 1))
  )
bi_thr <- fx$cfg$blood_cor_max
lab_bi <- filter(bi_df, class %in% c("removed: blood-tracking", "muscle marker"))
pal_bi <- c("removed: blood-tracking" = "#B2182B", "muscle marker" = "#2166AC", "ribosomal" = "#F1A340", "other" = "grey80")

blood_index_fig <- ggplot(bi_df, aes(blood_cor, ratio)) +
  annotate("rect", xmin = bi_thr, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#B2182B", alpha = 0.06) +
  geom_point(data = ~ filter(.x, class == "other"), colour = "grey82", size = 0.7, alpha = 0.4) +
  geom_point(data = ~ filter(.x, class != "other"), aes(colour = class), size = 1.9, alpha = 0.9) +
  geom_vline(xintercept = bi_thr, linetype = "dashed", colour = "#B2182B", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = lab_bi, aes(label = gene, colour = class), size = 2.1,
    max.overlaps = 16, seed = 42, min.segment.length = 0, segment.size = 0.2, show.legend = FALSE
  ) +
  scale_colour_manual(values = pal_bi, name = NULL) +
  labs(
    title = "Data-driven blood filter vs the erythrocyte-RNA ratio",
    subtitle = sprintf(
      "Blood-index correlation (x, cut at %.2f) removes red-cell without touching muscle; HPA's RNA ratio (y) would trap ribosomes",
      bi_thr
    ),
    x = "Spearman correlation with the per-sample blood index",
    y = expression(log[2] ~ "(erythrocyte / myonuclei RNA), HPA")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey35"),
    panel.grid.minor = element_blank(),
    legend.position = c(0.16, 0.16),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
  )

ggsave(file.path(report_dir, "filtering_blood_index.pdf"), blood_index_fig, width = 8, height = 6)
ggsave(file.path(report_dir, "filtering_blood_index.png"), blood_index_fig, width = 8, height = 6, dpi = 200)
cat("wrote", file.path(report_dir, "filtering_blood_index.png"), "\n")
