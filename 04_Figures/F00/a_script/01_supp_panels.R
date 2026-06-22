#!/usr/bin/env Rscript
# F00 — Pipeline QC Supplementary Figures for HRvLR
#
# Optional stitched QC composite built from canonical stage-local outputs:
#   page 1: normalization and preprocessing QC
#   page 2: imputation and DEP robustness QC

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

source("04_Figures/shared/style.R")

BASE <- "04_Figures/F00"
RPT <- file.path(BASE, "b_reports", "supp")
PNL <- file.path(RPT, "panels")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT, PNL, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

PW <- 89
PH <- 65

int_norm <- readRDS("02_Normalization/c_data/00_report_intermediates.rds")
norm_dal <- readRDS("02_Normalization/c_data/DAList_normalized.rds")
imp_dal <- readRDS("02_Normalization/imputation/c_data/DAList_imputed_missforest.rds")
da_summ <- read_csv("03_DEP/a_non_imputed/c_data/02_DA_summary.csv", show_col_types = FALSE)
outlier_sens <- if (file.exists("03_DEP/a_non_imputed/c_data/11_outlier_sensitivity.csv")) {
  read_csv("03_DEP/a_non_imputed/c_data/11_outlier_sensitivity.csv", show_col_types = FALSE)
} else {
  tibble()
}
imputation_sens <- if (file.exists("03_DEP/a_non_imputed/c_data/13_imputation_sensitivity.csv")) {
  read_csv("03_DEP/a_non_imputed/c_data/13_imputation_sensitivity.csv", show_col_types = FALSE)
} else {
  tibble()
}

contrast_order <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_HR", "Training_LR", "Acute_HR", "Acute_LR",
  "Training_Interaction", "Acute_Interaction"
)

save_png_pdf <- function(plot, stem, width_mm, height_mm) {
  save_panel(plot, file.path(PNL, stem), width_mm, height_mm)
}

compact_panel <- function(plot) {
  plot +
    labs(subtitle = NULL) +
    theme(
      plot.title = element_text(face = "bold", size = 7, margin = margin(b = 1)),
      legend.position = "bottom",
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5.5),
      legend.margin = margin(0, 0, 0, 0),
      plot.tag = element_text(face = "bold", size = 12)
    )
}

placeholder_plot <- function(tag, title, subtitle) {
  ggplot() +
    annotate("text",
      x = 0.5, y = 0.55, label = title,
      size = 4, fontface = "bold", color = "grey25"
    ) +
    annotate("text",
      x = 0.5, y = 0.42, label = subtitle,
      size = 3, color = "grey40"
    ) +
    labs(tag = tag) +
    theme_void() +
    theme(plot.tag = element_text(face = "bold", size = 12))
}

# Panel A: filter cascade ------------------------------------------------------

fcasc <- int_norm$filter_log %>%
  mutate(
    step = factor(step, levels = step),
    removed = ifelse(is.na(n_removed), 0L, n_removed)
  )

pA <- ggplot(fcasc, aes(step, n_after)) +
  geom_col(fill = "#2166AC", width = 0.7) +
  geom_text(aes(label = comma(n_after)), vjust = -1.5, size = 2.2, fontface = "bold") +
  geom_text(aes(label = ifelse(removed > 0, paste0("-", comma(removed)), "")),
    vjust = -0.2, size = 2.0, color = "#B2182B"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
  labs(
    x = NULL, y = "Proteins retained", tag = "A",
    title = "Protein filter cascade",
    subtitle = sprintf("%s -> %s proteins", comma(fcasc$n_after[1]), comma(tail(fcasc$n_after, 1)))
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 5.5))

write_csv(fcasc, file.path(DAT, "panel_A_filter_cascade.csv"))
save_png_pdf(pA, "SUPP_panel_A_filter_cascade", PW, PH)

# Panel B: normalization method ranking ---------------------------------------

norm_scores <- int_norm$norm_scores %>%
  arrange(composite) %>%
  mutate(method = factor(method, levels = rev(method)))

pB <- ggplot(norm_scores, aes(composite, method, fill = method == "cycloess")) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", composite)),
    hjust = -0.1,
    size = 2.0, fontface = "bold"
  ) +
  scale_fill_manual(values = c(`TRUE` = "#E41A1C", `FALSE` = "grey70"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = "Composite rank score", y = NULL, tag = "B",
    title = "Normalization score ranking",
    subtitle = "Lower score is better; cycloess selected"
  ) +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 6))

write_csv(norm_scores, file.path(DAT, "panel_B_normalization_scores.csv"))
save_png_pdf(pB, "SUPP_panel_B_normalization_scores", PW, PH)

# Panel C: pre-normalization PCA ----------------------------------------------

pca_pre <- int_norm$pca_pre$scores %>%
  mutate(
    Group_Time = factor(Group_Time, levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")),
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  )

pC <- ggplot(pca_pre, aes(PC1, PC2, color = Group_Time, shape = Timepoint, fill = Group_Time)) +
  stat_ellipse(aes(group = Group_Time),
    geom = "polygon", alpha = 0.10,
    level = 0.80, linewidth = 0.25, show.legend = FALSE
  ) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_color_manual(values = PCA_COLORS, name = NULL) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = SHAPE_TP, name = NULL) +
  labs(
    x = sprintf("PC1 (%.1f%%)", int_norm$pca_pre$var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", int_norm$pca_pre$var_exp[2]),
    tag = "C",
    title = "PCA before normalization"
  ) +
  FIG_THEME +
  theme(legend.position = "top")

write_csv(pca_pre, file.path(DAT, "panel_C_pca_pre.csv"))
save_png_pdf(pC, "SUPP_panel_C_pca_pre", PW, PH)

# Panel D: post-normalization PCA ---------------------------------------------

pca_post <- int_norm$pca_post$scores %>%
  mutate(
    Group_Time = factor(Group_Time, levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")),
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  )

pD <- ggplot(pca_post, aes(PC1, PC2, color = Group_Time, shape = Timepoint, fill = Group_Time)) +
  stat_ellipse(aes(group = Group_Time),
    geom = "polygon", alpha = 0.10,
    level = 0.80, linewidth = 0.25, show.legend = FALSE
  ) +
  geom_point(size = 1.5, alpha = 0.85) +
  scale_color_manual(values = PCA_COLORS, name = NULL) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = SHAPE_TP, name = NULL) +
  labs(
    x = sprintf("PC1 (%.1f%%)", int_norm$pca_post$var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", int_norm$pca_post$var_exp[2]),
    tag = "D",
    title = "PCA after normalization"
  ) +
  FIG_THEME +
  theme(legend.position = "top")

write_csv(pca_post, file.path(DAT, "panel_D_pca_post.csv"))
save_png_pdf(pD, "SUPP_panel_D_pca_post", PW, PH)

# Panel E: outlier consensus ---------------------------------------------------

od <- int_norm$outlier_diag %>%
  mutate(
    status = case_when(
      consensus_outlier ~ "Consensus outlier",
      n_flags > 0 ~ "Flagged",
      TRUE ~ "Clean"
    ),
    status = factor(status, levels = c("Clean", "Flagged", "Consensus outlier"))
  )

pE <- ggplot(od, aes(mahal_dist, n_flags, color = status, shape = Timepoint)) +
  geom_jitter(width = 0.04, height = 0.10, size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c(
    "Clean" = "grey55",
    "Flagged" = "#E6A100",
    "Consensus outlier" = "#B2182B"
  ), name = NULL) +
  scale_shape_manual(values = SHAPE_TP, guide = "none") +
  scale_y_continuous(breaks = 0:4) +
  labs(
    x = "Mahalanobis distance", y = "QC flags", tag = "E",
    title = "Outlier detection consensus",
    subtitle = sprintf("%d samples with >=1 flag", sum(od$n_flags > 0, na.rm = TRUE))
  ) +
  FIG_THEME +
  theme(legend.position = "top")

write_csv(od, file.path(DAT, "panel_E_outlier_consensus.csv"))
save_png_pdf(pE, "SUPP_panel_E_outlier_consensus", PW, PH)

# Panel F: per-sample missingness ---------------------------------------------

miss_bar <- int_norm$miss_bar_data %>%
  filter(status == "Missing") %>%
  mutate(
    Group_Time = factor(Group_Time, levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3"))
  )

pF <- ggplot(miss_bar, aes(reorder(Col_ID, -n), n, fill = Group_Time, alpha = is_outlier)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = PCA_COLORS, name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 1), guide = "none") +
  labs(
    x = NULL, y = "Missing proteins", tag = "F",
    title = "Per-sample missingness",
    subtitle = sprintf("%d analysis samples after QC", dplyr::n_distinct(pca_post$Col_ID))
  ) +
  FIG_THEME +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")

write_csv(miss_bar, file.path(DAT, "panel_F_sample_missingness.csv"))
save_png_pdf(pF, "SUPP_panel_F_sample_missingness", PW, PH)

# Panel G: MAR/MNAR classification — not available under minimal imputation
# Per-protein MAR/MNAR classification was dropped when stage 02 moved to the
# CvH-style method arms (imp4p / MsCoreUtils / missForest).

pG <- placeholder_plot(
  "G", "MAR/MNAR classification not available",
  "Per-protein classification dropped under minimal imputation"
)
save_png_pdf(pG, "SUPP_panel_G_missingness_scatter", PW, PH)

# Panel H: classification counts — not available under minimal imputation

pH <- placeholder_plot(
  "H", "Missingness class counts not available",
  "Per-protein classification dropped under minimal imputation"
)
save_png_pdf(pH, "SUPP_panel_H_missingness_counts", PW, PH)

# Panel I: observed vs imputed density ----------------------------------------

norm_mat <- as.matrix(norm_dal$data)
imp_mat <- as.matrix(imp_dal$data)[rownames(norm_mat), colnames(norm_mat)]
was_na <- is.na(norm_mat)
oob_error <- imp_dal$imputation$oob_error
obs_vals <- as.numeric(norm_mat[!was_na])
imp_vals <- as.numeric(imp_mat[was_na])
dens_df <- bind_rows(
  tibble(value = obs_vals, type = "Observed"),
  tibble(value = imp_vals, type = "Imputed")
) %>%
  mutate(type = factor(type, levels = c("Observed", "Imputed")))

pI <- ggplot(dens_df, aes(value, fill = type, color = type)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  scale_fill_manual(values = c("Observed" = "#377EB8", "Imputed" = "#E41A1C"), name = NULL) +
  scale_color_manual(values = c("Observed" = "#377EB8", "Imputed" = "#E41A1C"), name = NULL) +
  annotate("text",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
    label = sprintf("OOB = %.3f", oob_error),
    size = 2.2, fontface = "italic", color = "grey30"
  ) +
  labs(
    x = "log2 intensity", y = "Density", tag = "I",
    title = "Observed vs imputed intensity"
  ) +
  FIG_THEME +
  theme(legend.position = "top")

write_csv(
  dens_df %>% group_by(type) %>% summarise(n = n(), mean = mean(value), median = median(value), sd = sd(value), .groups = "drop"),
  file.path(DAT, "panel_I_imputation_density_summary.csv")
)
save_png_pdf(pI, "SUPP_panel_I_imputation_density", PW, PH)

# Panel J: MNAR shift — not available under minimal imputation
# Mechanism-agnostic missForest imputes all missing values with one model and
# produces no MNAR-specific shift audit.

pJ <- placeholder_plot(
  "J", "MNAR shift audit not available",
  "Mechanism-agnostic missForest produces no MNAR audit"
)
save_png_pdf(pJ, "SUPP_panel_J_mnar_shift", PW, PH)

# Panel K: DEP counts heatmap --------------------------------------------------

dep_counts <- da_summ %>%
  filter(type %in% c("up", "down")) %>%
  group_by(contrast) %>%
  summarise(
    `p<0.10` = sum(sig.PVal, na.rm = TRUE),
    `FDR<0.10` = sum(sig.FDR.10, na.rm = TRUE),
    `FDR<0.05` = sum(sig.FDR.05, na.rm = TRUE),
    `Pi<0.05` = sum(sig.Pi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-contrast, names_to = "threshold", values_to = "n") %>%
  mutate(
    contrast = factor(contrast, levels = intersect(contrast_order, unique(contrast))),
    threshold = factor(threshold, levels = c("p<0.10", "FDR<0.10", "FDR<0.05", "Pi<0.05"))
  )

pK <- ggplot(dep_counts, aes(threshold, contrast, fill = n)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = comma(n), color = n > max(n, na.rm = TRUE) * 0.5),
    size = 2.3, fontface = "bold"
  ) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  scale_fill_gradient(low = "#DEEBF7", high = "#08519C", name = "DEPs") +
  labs(
    x = NULL, y = NULL, tag = "K",
    title = "DEP counts by contrast and threshold"
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1, size = 6),
    legend.position = "right"
  )

write_csv(dep_counts, file.path(DAT, "panel_K_dep_counts.csv"))
save_png_pdf(pK, "SUPP_panel_K_dep_heatmap", PW, PH)

# Panel L: stage-03 robustness -------------------------------------------------

robust_long <- bind_rows(
  if (nrow(outlier_sens) > 0) {
    outlier_sens %>%
      transmute(
        contrast = Contrast,
        check = "Outlier reduction",
        rho = Spearman_rho,
        delta_pi = Pi_reduced - Pi_full
      )
  } else {
    tibble()
  },
  if (nrow(imputation_sens) > 0) {
    imputation_sens %>%
      transmute(
        contrast = contrast,
        check = "Imputation sensitivity",
        rho = spearman_rho,
        delta_pi = NA_real_
      )
  } else {
    tibble()
  }
) %>%
  mutate(
    contrast = factor(contrast, levels = intersect(contrast_order, unique(contrast))),
    check = factor(check, levels = c("Outlier reduction", "Imputation sensitivity"))
  )

if (nrow(robust_long) > 0) {
  pi_labels <- robust_long %>%
    filter(check == "Outlier reduction", !is.na(delta_pi)) %>%
    mutate(label = sprintf("dPi=%+d", as.integer(round(delta_pi))))

  pL <- ggplot(robust_long, aes(contrast, rho, color = check, group = check)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    geom_point(size = 2) +
    geom_text(
      data = pi_labels, aes(label = label), inherit.aes = FALSE,
      x = pi_labels$contrast, y = pi_labels$rho + 0.0025,
      size = 2.1, color = "grey25", fontface = "bold"
    ) +
    scale_color_manual(values = c("Outlier reduction" = "#2166AC", "Imputation sensitivity" = "#B2182B"), name = NULL) +
    scale_y_continuous(limits = c(0.94, 1.005), breaks = c(0.95, 0.97, 0.99, 1.00)) +
    labs(
      x = NULL, y = "Spearman rho", tag = "L",
      title = "Stage-03 robustness checks",
      subtitle = "Rank concordance after stage perturbation"
    ) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
      legend.position = "top"
    )
  write_csv(robust_long, file.path(DAT, "panel_L_robustness_summary.csv"))
} else {
  pL <- placeholder_plot("L", "Robustness summaries unavailable", "Sensitivity tables were not found")
}

save_png_pdf(pL, "SUPP_panel_L_robustness", PW, PH)

# Composite assembly -----------------------------------------------------------

page1 <- (compact_panel(pA) | compact_panel(pB)) /
  (compact_panel(pC) | compact_panel(pD)) /
  (compact_panel(pE) | compact_panel(pF)) +
  plot_layout(heights = c(0.9, 1, 1)) +
  plot_annotation(
    title = "Pipeline QC - Pre-Processing",
    subtitle = sprintf(
      "HPA filter, contaminant removal, cycloess normalization, and consensus QC | %s -> %s proteins",
      comma(int_norm$n_raw), comma(tail(fcasc$n_after, 1))
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "italic", size = 8, color = "grey30")
    )
  )

page2 <- (compact_panel(pG) | compact_panel(pH)) /
  (compact_panel(pI) | compact_panel(pJ)) /
  (compact_panel(pK) | compact_panel(pL)) +
  plot_layout(heights = c(0.9, 1, 1)) +
  plot_annotation(
    title = "Pipeline QC - Imputation and DEP",
    subtitle = sprintf(
      "MAR/MNAR classification, missForest (OOB = %.3f), and DEP robustness summaries",
      oob_error
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "italic", size = 8, color = "grey30")
    )
  )

ggsave(file.path(RPT, "SUPP_F00_normalization.pdf"), page1,
  width = 178, height = 205, units = "mm", device = PDF_DEVICE
)
ggsave(file.path(RPT, "SUPP_F00_normalization.png"), page1,
  width = 178, height = 205, units = "mm", dpi = 300
)

ggsave(file.path(RPT, "SUPP_F00_imputation_dep.pdf"), page2,
  width = 178, height = 205, units = "mm", device = PDF_DEVICE
)
ggsave(file.path(RPT, "SUPP_F00_imputation_dep.png"), page2,
  width = 178, height = 205, units = "mm", dpi = 300
)

message("F00 supplementary QC composites done")
