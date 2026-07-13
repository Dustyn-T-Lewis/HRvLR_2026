#!/usr/bin/env Rscript
# Filtering QC figure. Reports what was removed and what contamination remained, rather than
# deleting the evidence and asserting the samples were clean.

pacman::p_load(here, readr, dplyr, tidyr, stringr, ggplot2, forcats, patchwork, scales)

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
