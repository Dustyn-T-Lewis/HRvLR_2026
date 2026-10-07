# Supplement to Figure 1: data and analysis quality. The top row asks whether the samples held up
# (filtering, detection, intensity before and after normalisation, dynamic range, missingness);
# the bottom row whether the analysis behaves (sample correlation, blood, blood adjustment, p-value
# calibration, power). Reads saved outputs only; the summaries here (ordering, within- and
# between-subject correlation, histograms) describe them and test nothing.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source(here("R", "helpers.R"))
source(here("R", "figure_theme.R"))
source(here("R", "figure_panels.R"))

stage <- here("01_Trunk", "Figures", "S01_QC")
figure_dir <- file.path(stage, "b_reports")
out <- file.path(stage, "c_data")
walk(c(figure_dir, out), dir.create, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  filter = "01_Trunk/01_Protein_Filtering/c_data/protein_filtering.xlsx",
  normalize = "01_Trunk/02_Normalisation/c_data/normalisation.xlsx",
  quantify = "01_Trunk/03_Quantification/c_data/quantification.xlsx",
  contrasts = "02_HRvLR_groups/02_Protein_Contrasts/c_data/protein_contrasts.xlsx"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Missing inputs: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
read_sheet <- \(input, sheet) readxl::read_excel(paths[[input]], sheet)
contrast_names <- c(
  Training_HR = "Training HR", Training_LR = "Training LR", Acute_HR = "Acute HR",
  Acute_LR = "Acute LR", Baseline_HRvLR = "HR - LR T1", Trained_HRvLR = "HR - LR T2",
  Acute_HRvLR = "HR - LR T3", Training_Interaction = "Int. training",
  Acute_Interaction = "Int. acute"
)
contrast_family <- c(
  Training_HR = "Training", Training_LR = "Training", Acute_HR = "Acute", Acute_LR = "Acute",
  Baseline_HRvLR = "HR - LR", Trained_HRvLR = "HR - LR", Acute_HRvLR = "HR - LR",
  Training_Interaction = "Interaction", Acute_Interaction = "Interaction"
)
time_axis <- scale_colour_manual(values = time_colours, name = NULL)

# a: proteins left after each filtering step.
filter_log <- read_sheet("filter", "filter_log") |>
  mutate(step = factor(step, rev(step)))
panel_filter <- ggplot(filter_log, aes(n_after, step)) +
  geom_col(fill = "grey70", width = 0.7) +
  geom_text(aes(label = n_after), hjust = 1.1, size = 1.8, colour = "white") +
  labs(
    x = "proteins remaining", y = NULL, title = "Filtering",
    subtitle = sprintf("%d to %d proteins", max(filter_log$n_after), min(filter_log$n_after))
  ) +
  figure_theme +
  theme(panel.grid.major.y = element_blank())

# b: share of proteins missing in each of the 48 runs. The three runs not carried forward are the
# ones absent from the 45-sample blood table that 03_Quantification writes.
kept_samples <- read_sheet("quantify", "blood")$sample_id
detection <- read_sheet("filter", "outlier_diagnostics") |>
  mutate(
    kept = sample_id %in% kept_samples,
    sample_id = factor(sample_id, sample_id[order(arm, timepoint, pct_missing)])
  )
panel_detection <- ggplot(detection, aes(sample_id, pct_missing, colour = timepoint)) +
  geom_segment(aes(xend = sample_id, y = 0, yend = pct_missing), linewidth = 0.3) +
  geom_point(aes(shape = kept), size = 1.2) +
  facet_grid(~arm, scales = "free_x", space = "free_x") +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 4), labels = c("excluded", "kept"), name = NULL
  ) +
  time_axis +
  labs(
    x = "run", y = "% of proteins missing", title = "Detection per run",
    subtitle = sprintf("%d runs; %d excluded as outliers", nrow(detection), sum(!detection$kept))
  ) +
  figure_theme +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())

# c: each sample's intensity quantiles before and after cyclic loess.
quantiles <- read_sheet("normalize", "intensity_quantiles") |>
  mutate(
    stage = factor(stage, c("before", "after"), c("before normalisation", "after")),
    sample_id = factor(sample_id, unique(sample_id[order(arm, timepoint)]))
  )
panel_intensity <- ggplot(quantiles, aes(sample_id, colour = timepoint)) +
  geom_linerange(aes(ymin = q05, ymax = q95), linewidth = 0.25) +
  geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.9) +
  geom_point(aes(y = q50), colour = "white", size = 0.4) +
  facet_wrap(~stage, ncol = 1) +
  time_axis +
  labs(
    x = "sample", y = "log2 intensity", title = "Intensity per sample",
    subtitle = "median, 25-75% and 5-95% of proteins"
  ) +
  figure_theme +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())

# d: dynamic range, with abundant muscle proteins named.
protein_table <- read_sheet("quantify", "proteins") |>
  arrange(desc(mean_complete)) |>
  mutate(rank = row_number())
landmarks <- c("MYH7", "MYH2", "ACTA1", "TTN", "DES", "CKM", "MB", "ATP5F1A", "NDUFA4", "PSMA1")
named <- filter(protein_table, gene %in% landmarks)
panel_range <- ggplot(protein_table, aes(rank, mean_complete)) +
  geom_point(size = 0.3, colour = "grey60") +
  geom_point(data = named, size = 1, colour = "black") +
  ggrepel::geom_text_repel(
    aes(label = gene), named,
    size = 1.8, min.segment.length = 0, segment.size = 0.2, seed = 1
  ) +
  labs(
    x = "rank by mean abundance", y = "mean log2 abundance", title = "Dynamic range",
    subtitle = sprintf(
      "%d proteins over %.1f log2 units",
      nrow(protein_table), diff(range(protein_table$mean_complete))
    )
  ) +
  figure_theme

# e: share missing against mean observed abundance, with limpa's detection probability curve.
curve <- read_sheet("quantify", "curve")
dpc <- tibble(mean_observed = seq(
  min(protein_table$mean_observed), max(protein_table$mean_observed),
  length.out = 200
)) |>
  mutate(missing = 1 - plogis(curve$estimate[1] + curve$estimate[2] * mean_observed))
panel_missing <- ggplot(protein_table, aes(mean_observed, 1 - prop_observed)) +
  geom_point(size = 0.3, alpha = 0.4, colour = "grey40") +
  geom_line(aes(y = missing), dpc, colour = contrast_tone("line")[["HR - LR"]], linewidth = 0.6) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    x = "mean log2 abundance where observed", y = "values missing",
    title = "Missing because faint",
    subtitle = sprintf(
      "%.1f%% missing; line: limpa detection curve, slope %.2f",
      curve$pct_missing[1], curve$estimate[2]
    )
  ) +
  figure_theme

# f: correlation between the 45 samples, each protein centred, ordered by clustering.
correlation <- read_sheet("quantify", "sample_correlation")
r <- as.matrix(column_to_rownames(correlation, "sample_id"))
ordering <- rownames(r)[hclust(as.dist(1 - r), method = "average")$order]
subject_of <- setNames(read_sheet("quantify", "blood")$subject, kept_samples)
pairs <- as_tibble(r, rownames = "a") |>
  pivot_longer(-a, names_to = "b", values_to = "r") |>
  filter(a < b) |>
  mutate(same = subject_of[a] == subject_of[b])
within_between <- summarise(pairs, r = median(r), .by = same)
cells <- as_tibble(r, rownames = "a") |>
  pivot_longer(-a, names_to = "b", values_to = "r") |>
  # The diagonal (r = 1) is left blank so it does not set the colour range.
  mutate(r = if_else(a == b, NA_real_, r), a = factor(a, ordering), b = factor(b, ordering))
panel_correlation <- ggplot(cells, aes(a, b, fill = r)) +
  geom_tile() +
  scale_fill_distiller(
    palette = "RdBu", limits = c(-0.6, 0.6), breaks = c(-0.5, 0, 0.5), name = "r",
    oob = scales::squish, na.value = "white"
  ) +
  labs(
    x = NULL, y = NULL, title = "Sample correlation",
    subtitle = sprintf(
      "median r: same person %.2f, different people %.2f",
      within_between$r[within_between$same], within_between$r[!within_between$same]
    )
  ) +
  figure_theme +
  theme(axis.text = element_blank(), panel.grid = element_blank())

# g: blood index per sample, by arm and time.
blood <- read_sheet("quantify", "blood") |>
  mutate(arm = sub("_.*", "", group), timepoint = sub(".*_", "", group))
blood_test <- filter(read_sheet("filter", "blood_by_arm"), scale == "log2_index")
panel_blood <- ggplot(blood, aes(timepoint, blood_rel, colour = timepoint, shape = arm)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, seed = 1), size = 1.1) +
  scale_shape_manual(values = arm_shapes, name = NULL) +
  time_axis +
  labs(
    x = NULL, y = "blood index (log2, vs sample median)", title = "Blood in each biopsy",
    subtitle = sprintf(
      "LR - HR at T3 %.2f, permutation p = %.2g", blood_test$lr_minus_hr_at_t3,
      blood_test$permutation_p
    )
  ) +
  figure_theme

# h: each protein's log2 FC with and without the blood covariate, every contrast pooled.
dep <- read_sheet("contrasts", "DEP_matrix")
unadjusted <- read_sheet("contrasts", "unadjusted")
blood_effect <- names(contrast_family) |>
  map(\(cn) {
    tibble(
      contrast = cn, family = contrast_family[[cn]],
      adjusted = dep[[paste0(cn, "_logFC")]],
      unadjusted = unadjusted[[paste0(cn, "_logFC")]][match(dep$uniprot_id, unadjusted$uniprot_id)]
    )
  }) |>
  list_rbind()
blood_r <- cor(blood_effect$adjusted, blood_effect$unadjusted)
panel_blood_fit <- ggplot(blood_effect, aes(unadjusted, adjusted, colour = family)) +
  geom_abline(linetype = "dashed", linewidth = 0.25) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_colour_manual(values = contrast_tone("line"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(
    x = "log2 FC without blood covariate", y = "log2 FC with blood covariate",
    title = "Blood adjustment", subtitle = sprintf("nine contrasts pooled, r = %.3f", blood_r)
  ) +
  figure_theme

# i: p-value histograms; the dashed line is a flat null.
calibration <- read_sheet("contrasts", "null_calibration")
p_values <- names(contrast_family) |>
  map(\(cn) tibble(contrast = cn, family = contrast_family[[cn]], p = dep[[paste0(cn, "_p")]])) |>
  list_rbind() |>
  mutate(contrast = factor(contrast, names(contrast_family)))
panel_p <- ggplot(p_values, aes(p, fill = family)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  geom_hline(yintercept = nrow(dep) / 20, linetype = "dashed", linewidth = 0.25) +
  geom_text(
    aes(0.5, Inf, label = sprintf("median |t| %.2f", median_abs_t)),
    mutate(calibration, contrast = factor(contrast, names(contrast_family))),
    inherit.aes = FALSE, vjust = 1.3, size = 1.6
  ) +
  facet_wrap(~contrast, ncol = 3, labeller = as_labeller(contrast_names)) +
  scale_fill_manual(values = contrast_tone("fill"), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(
    x = "p", y = "proteins", title = "p-value calibration",
    subtitle = "dashed: flat null; median |t| is 0.67 under the null"
  ) +
  figure_theme +
  theme(
    strip.text = element_text(size = 5), panel.grid.minor = element_blank(),
    panel.spacing.x = unit(2.5, "mm")
  )

# j: observed effects in the primary contrast against the effect it can detect at 80% power.
detectable <- read_sheet("contrasts", "detectable_effect")
observed <- tibble(effect = abs(dep$Training_Interaction_logFC))
panel_power <- ggplot(observed, aes(effect)) +
  geom_histogram(binwidth = 0.1, boundary = 0, fill = contrast_tone("fill")[["Interaction"]]) +
  geom_vline(
    aes(xintercept = log2_fc, linetype = factor(quantile)), detectable,
    linewidth = 0.4, colour = contrast_tone("line")[["Interaction"]]
  ) +
  scale_linetype_manual(
    values = c(`0.5` = "solid", `0.9` = "dashed"),
    labels = c("median protein", "90th percentile"), name = "detectable at 80% power"
  ) +
  labs(
    x = "|log2 FC|, training interaction", y = "proteins", title = "Power",
    subtitle = sprintf(
      "%.0f%% of effects below the median detectable %.2f",
      100 * mean(observed$effect < detectable$log2_fc[detectable$quantile == 0.5]),
      detectable$log2_fc[detectable$quantile == 0.5]
    )
  ) +
  figure_theme +
  theme(legend.position = "bottom", legend.title.position = "top")

# Every plotted number traces to its source table before anything is written.
stopifnot(
  min(filter_log$n_after) == nrow(protein_table),
  nrow(detection) == 48, sum(!detection$kept) == 3, length(kept_samples) == 45,
  nrow(quantiles) == 2 * 45, nrow(r) == 45, isSymmetric(r),
  all(!is.na(blood_effect$unadjusted)), nrow(p_values) == 9 * nrow(dep),
  nrow(calibration) == 9, nrow(detectable) == 2
)

figure <- wrap_plots(
  panel_filter, panel_detection, panel_intensity, panel_range, panel_missing,
  panel_correlation, panel_blood, panel_blood_fit, panel_p, panel_power,
  nrow = 2
) +
  plot_annotation(tag_levels = "a") &
  figure_tags &
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"))
save_figure(figure, figure_dir, "s04_qc", width = 300, height = 190)

sheets <- list(
  filter_log = select(filter_log, step, n_removed, n_after),
  detection = select(detection, sample_id, subject, arm, timepoint, pct_missing, kept),
  intensity_quantiles = select(quantiles, -arm, -timepoint),
  dynamic_range = select(protein_table, uniprot_id, gene, rank, mean_complete),
  missingness = select(protein_table, uniprot_id, gene, mean_observed, prop_observed),
  correlation_summary = within_between,
  blood = select(blood, sample_id, subject, arm, timepoint, blood_rel),
  blood_adjustment = tibble(r = blood_r, n_points = nrow(blood_effect)),
  calibration = calibration,
  power = detectable
)
holds <- c(
  "Panel a: proteins removed and remaining at each filtering step.",
  "Panel b: share of proteins missing per run, and whether the run was kept.",
  "Panel c: per-sample log2 intensity quantiles before and after cyclic loess.",
  "Panel d: proteins ranked by mean log2 abundance over the complete matrix.",
  "Panel e: per protein, mean observed abundance and share of values observed.",
  "Panel f: median sample correlation within and between people.",
  "Panel g: blood index per sample.",
  "Panel h: correlation of log2 FC with and without the blood covariate, nine contrasts pooled.",
  "Panel i: per contrast, propTrueNull, median |t| and share of p below 0.2.",
  "Panel j: log2 FC the training interaction detects at 80% power."
)
write_workbook(sheets, holds, file.path(out, "s04_qc.xlsx"), inputs, paths)
message("wrote s04_qc.pdf, s04_qc.png and s04_qc.xlsx")
