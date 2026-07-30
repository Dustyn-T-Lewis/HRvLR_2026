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
      paste(
        "Blood rises at T3 (hyperaemia) and leukocytes rise with it",
        "(post-exercise infiltration). The rise is NOT equal in the two arms: on"
      ),
      paste(
        "the log2 haemoglobin index, arm x T3 gives b = -1.21, p = 0.032, and",
        "p = 0.017 by subject-label permutation, LR rising 2.14 against HR's 0.57."
      ),
      paste(
        "On this percent-of-signal scale the same term is p = 0.26 - a bounded",
        "proportion with far larger relative residual variance, so it fails to"
      ),
      paste(
        "detect what the index measures. The confound does not cancel in the",
        "interaction. Keratin and adipose are negligible throughout."
      ),
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

# Supplement: the curated blood list against both lines of evidence. Each protein's correlation
# with the per-sample blood index (x) against HPA's erythrocyte/myonuclei RNA ratio (y), which on
# its own would trap the ribosomal proteins and miss the enucleate red-cell membrane skeleton.
muscle_markers <- c(
  "MB", "CKM", "ACTA1", "MYH1", "MYH2", "MYH7", "ALDOA", "CASQ1", "PYGM",
  "ATP2A1", "DES", "TNNT3", "TNNI2", "MYBPC1", "TNNC2", "MYL1"
)
# The panel now checks the curated blood list against the in-sample evidence
# instead of drawing the threshold that used to make the call. A curated protein
# sitting at blood_cor near zero would mean the list is wrong; a kept protein
# sitting high means covariation alone was never enough to remove it.
blood_list <- readr::read_csv(
  here("00_input", "blood_contaminants.csv"),
  show_col_types = FALSE
)
bi_df <- fx$protein_calls |>
  filter(!is.na(blood_cor)) |>
  mutate(
    class = case_when(
      gene %in% blood_list$gene ~ "removed: curated blood",
      gene %in% muscle_markers ~ "muscle marker",
      str_detect(gene, "^RP[LS]") ~ "ribosomal",
      TRUE ~ "other"
    ),
    ratio = log2((ery + 1) / (myo + 1))
  )

# 99.9th percentile of a 300x permutation of the index, restricted to proteins
# with enough observations to be callable (00_blood_cor_null.R). Drawn as
# context, not as a gate: nothing is removed for crossing it.
bi_thr <- 0.59
lab_bi <- filter(bi_df, class %in% c("removed: curated blood", "muscle marker"))
pal_bi <- c("removed: curated blood" = "#B2182B", "muscle marker" = "#2166AC", "ribosomal" = "#F1A340", "other" = "grey80")

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
    title = "The curated blood list against the in-sample blood index",
    subtitle = sprintf(
      "Curated blood proteins track the haemoglobin index; muscle markers do not. Dashed line is the permutation null (%.2f), context only \u2014 removal is by protein identity",
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
