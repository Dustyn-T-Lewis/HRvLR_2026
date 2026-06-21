#!/usr/bin/env Rscript
# 02_imputation_reports.R — Missingness and imputation quality reports
#
# HRvLR 2x3 design: Responder (HR/LR) x Timepoint (T1/T2/T3)
#
# Reads:
#   c_data/00_report_intermediates.rds
#   c_data/benchmark/04_composite_ranking.csv  (preferred, composite ranking)
#   c_data/benchmark/03_benchmark_summary.csv  (fallback, old summary)
#
# Writes:
#   b_reports/01_missingness_report.pdf  (1 page, 16x10 in)
#   b_reports/02_imputation_report.pdf   (2 pages, 14x10 in)

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(readr)
library(scales)

setwd(rprojroot::find_rstudio_root_file())

REPORT_DIR <- "02_Imputation/b_reports"
DATA_DIR   <- "02_Imputation/c_data"
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)

THM    <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank())

# --- Load intermediates ------------------------------------------------------
rpt <- readRDS(file.path(DATA_DIR, "00_report_intermediates.rds"))

mat            <- rpt$mat
mat_imp        <- rpt$mat_imp
was_na         <- rpt$was_na
meta           <- rpt$meta
miss_class     <- rpt$miss_class
miss_by_group  <- rpt$miss_by_group
prot_pct       <- rpt$prot_pct
pct_miss       <- rpt$pct_miss
mnar_genes     <- rpt$mnar_genes
mnar_audit     <- rpt$mnar_audit
best           <- rpt$best
n_mar_prots    <- rpt$n_mar_prots
n_mnar_prots   <- rpt$n_mnar_prots
mar_miss_vals  <- rpt$mar_miss_vals
mnar_miss_vals <- rpt$mnar_miss_vals
total_miss_vals <- rpt$total_miss_vals
classification_method <- rpt$classification_method
PAL_GT         <- rpt$PAL_GT
PAL_MAR        <- rpt$PAL_MAR
PAL_CLASS      <- rpt$PAL_CLASS

n_prot       <- nrow(mat)
n_samp       <- ncol(mat)
n_total_miss <- total_miss_vals
n_unreliable <- sum(!mnar_audit$imputation_reliable)

# =============================================================================
# Report 1: Missingness (1 page, 16x10)
# =============================================================================

# Panel A — Per-protein missingness histogram
pA1 <- miss_class |>
  filter(classification != "Complete") |>
  ggplot(aes(x = pct_miss, fill = classification)) +
  geom_histogram(binwidth = 5, boundary = 0, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = PAL_CLASS, name = NULL) +
  scale_x_continuous(labels = label_percent(scale = 1)) +
  labs(title = "A. Per-protein missingness", x = "% missing", y = "Count") +
  THM + theme(legend.position = "top")

# Panel B — Classification summary stacked bar
class_counts <- miss_class |>
  count(classification) |>
  mutate(classification = factor(classification, levels = c("Complete", "MAR", "MNAR")))

pB1 <- class_counts |>
  ggplot(aes(x = 1, y = n, fill = classification)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = paste0(classification, "\n", n)),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = PAL_CLASS) +
  coord_flip() +
  labs(title = "B. Classification summary", x = NULL, y = "Proteins") +
  THM + theme(legend.position = "none", axis.text.y = element_blank())

# Panel C — Scatter: mean intensity vs % missing
pC1 <- miss_class |>
  filter(classification != "Complete") |>
  ggplot(aes(x = mean_intensity, y = pct_miss, color = classification)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = PAL_CLASS, name = NULL) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(title = "C. Intensity vs missingness",
       x = "Mean log2 intensity", y = "% missing") +
  THM + theme(legend.position = "top")

# Panel D — Per-sample missingness bar
sample_miss <- data.frame(
  sample = colnames(mat),
  n_miss = colSums(is.na(mat))
) |>
  left_join(meta, by = c("sample" = "Col_ID")) |>
  mutate(Group_Time = factor(Group_Time, levels = names(PAL_GT)))

pD1 <- sample_miss |>
  ggplot(aes(x = reorder(sample, n_miss), y = n_miss, fill = Group_Time)) +
  geom_col() +
  scale_fill_manual(values = PAL_GT, name = NULL) +
  coord_flip() +
  labs(title = "D. Per-sample missingness", x = NULL, y = "Missing proteins") +
  THM + theme(legend.position = "top",
              axis.text.y = element_text(size = 6))

# Page 2 of missingness: per-group heatmap (top 50)
# NB: prot_pct + miss_by_group are saved in original (pre-gene_order) row order;
# rownames(mat) in rpt is gene-order sorted. Use miss_by_group's rownames here.
top_idx <- order(prot_pct, decreasing = TRUE)[1:min(50, sum(prot_pct > 0))]
p_heat <- as_tibble(miss_by_group[top_idx, ], rownames = "gene") |>
  pivot_longer(-gene, names_to = "Group_Time", values_to = "pct") |>
  mutate(gene = factor(gene, levels = rev(rownames(miss_by_group)[top_idx]))) |>
  ggplot(aes(Group_Time, gene, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "white", mid = "#FDDBC7", high = "#B2182B",
                       midpoint = 50, name = "% Miss") +
  labs(x = NULL, y = NULL, title = "E. Per-group missingness (top 50)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1))

class_method <- classification_method %||% "k-means"
miss_sub <- sprintf("%s proteins \u00d7 %s samples | %.1f%% missing | %s classification",
                    comma(n_prot), n_samp, pct_miss, class_method)

fig1 <- (pA1 | pB1) / (pC1 | pD1) +
  plot_annotation(
    title = "Missingness Report",
    subtitle = miss_sub,
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

pdf(file.path(REPORT_DIR, "01_missingness_report.pdf"), width = 16, height = 10)
print(fig1)
print(
  p_heat +
    plot_annotation(title = "Per-Group Missingness Detail"))
dev.off()
cat("Wrote:", file.path(REPORT_DIR, "01_missingness_report.pdf"), "\n")

# =============================================================================
# Report 2: Imputation quality (2 pages, 14x10)
# =============================================================================

# --- Page 1: Benchmark ranking (if available) or OOB summary ----------------

# Read OOB from summary file
oob_line <- grep("oob_error", readLines(file.path(DATA_DIR, "09_imputation_summary.txt")),
                 value = TRUE)
oob_val <- as.numeric(sub(".*=\\s*(?:c\\(NRMSE = )?([0-9.]+)\\)?", "\\1", oob_line, perl = TRUE))

# Try loading benchmark ranking — prefer new composite, fall back to old summary
composite_path <- file.path(DATA_DIR, "benchmark", "04_composite_ranking.csv")
bench_path <- file.path(DATA_DIR, "benchmark", "03_benchmark_summary.csv")
has_composite <- file.exists(composite_path)
has_bench <- has_composite || file.exists(bench_path)

if (has_composite) {
  comp_df <- read_csv(composite_path, show_col_types = FALSE)

  pA2 <- comp_df |>
    mutate(is_mf = grepl("^missForest|^RF", method),
           method = factor(method, levels = rev(method[order(composite)]))) |>
    ggplot(aes(x = composite, y = method, alpha = is_mf)) +
    geom_point(size = 3.5, color = "#2166AC") +
    scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.4), guide = "none") +
    coord_flip() +
    labs(x = "Composite score (higher = better)", y = NULL,
         title = "A. Composite benchmark ranking",
         subtitle = sprintf("Best: %s (%.3f) | %d methods",
                            comp_df$method[1], comp_df$composite[1],
                            nrow(comp_df))) + THM

  samp_means <- data.frame(
    Col_ID = colnames(mat),
    Pre  = colMeans(mat, na.rm = TRUE),
    Post = colMeans(mat_imp)) |>
    left_join(meta |> select(Col_ID, Group_Time), by = "Col_ID")

  pB2 <- ggplot(samp_means, aes(Pre, Post, color = Group_Time)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = PAL_GT) +
    labs(x = "Pre-imputation sample mean", y = "Post-imputation sample mean",
         title = "B. Sample-level mean shift",
         subtitle = sprintf("OOB error = %.4f", oob_val)) +
    THM + theme(legend.position = "bottom")

  page1 <- (pA2 | pB2) +
    plot_annotation(
      title = "Imputation Benchmark",
      subtitle = sprintf("%d methods | composite ranking | OOB = %.4f",
                         nrow(comp_df), oob_val),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

} else if (file.exists(bench_path)) {
  bench_sum <- read_csv(bench_path, show_col_types = FALSE)

  visible  <- bench_sum |> filter(median_nrmse < 5)
  excluded <- setdiff(bench_sum$method, visible$method)
  excl_note <- if (length(excluded) > 0)
    paste0(" | Off-scale: ", paste(excluded, collapse = ", ")) else ""

  pA2 <- visible |>
    ggplot(aes(x = reorder(method, -median_nrmse), y = median_nrmse)) +
    geom_point(size = 3.5, color = "#2166AC") +
    geom_errorbar(aes(ymin = median_nrmse - sd_nrmse,
                      ymax = median_nrmse + sd_nrmse),
                  width = 0.3, linewidth = 0.5, color = "#2166AC") +
    coord_flip() +
    labs(x = NULL, y = "NRMSE (lower = better)",
         title = "A. Benchmark NRMSE ranking",
         subtitle = sprintf("Best: %s (%.4f)%s",
                            bench_sum$method[which.min(bench_sum$median_nrmse)],
                            min(bench_sum$median_nrmse), excl_note)) + THM

  samp_means <- data.frame(
    Col_ID = colnames(mat),
    Pre  = colMeans(mat, na.rm = TRUE),
    Post = colMeans(mat_imp)) |>
    left_join(meta |> select(Col_ID, Group_Time), by = "Col_ID")

  pB2 <- ggplot(samp_means, aes(Pre, Post, color = Group_Time)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = PAL_GT) +
    labs(x = "Pre-imputation sample mean", y = "Post-imputation sample mean",
         title = "B. Sample-level mean shift",
         subtitle = sprintf("OOB error = %.4f", oob_val)) +
    THM + theme(legend.position = "bottom")

  page1 <- (pA2 | pB2) +
    plot_annotation(
      title = "Imputation Benchmark",
      subtitle = sprintf("missForest applied | OOB = %.4f", oob_val),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
} else {
  # No benchmark available — show OOB + sample shift only
  samp_means <- data.frame(
    Col_ID = colnames(mat),
    Pre  = colMeans(mat, na.rm = TRUE),
    Post = colMeans(mat_imp)) |>
    left_join(meta |> select(Col_ID, Group_Time), by = "Col_ID")

  pA2 <- ggplot(samp_means, aes(Pre, Post, color = Group_Time)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = PAL_GT) +
    labs(x = "Pre-imputation sample mean", y = "Post-imputation sample mean",
         title = "A. Sample-level mean shift",
         subtitle = sprintf("missForest OOB error = %.4f", oob_val)) +
    THM + theme(legend.position = "bottom")

  page1 <- pA2 +
    plot_annotation(
      title = "Imputation Summary",
      subtitle = sprintf("missForest | OOB = %.4f | No benchmark available", oob_val),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
}

# --- Page 2: Density overlay + MNAR audit + Cohen's d -----------------------

# Panel A — Density: observed vs imputed
obs_vals <- as.numeric(mat[!was_na])
imp_vals <- as.numeric(mat_imp[was_na])

dens_df <- bind_rows(
  data.frame(value = obs_vals, type = "Observed"),
  data.frame(value = imp_vals, type = "Imputed")
)

pA3 <- dens_df |>
  ggplot(aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c(Observed = "#377EB8", Imputed = "#E41A1C"), name = NULL) +
  scale_color_manual(values = c(Observed = "#377EB8", Imputed = "#E41A1C"), name = NULL) +
  labs(title = "A. Observed vs imputed value distributions",
       x = "log2 intensity", y = "Density") +
  THM + theme(legend.position = "top")

# Panel B — MNAR audit scatter: pre vs post mean
pB3 <- mnar_audit |>
  ggplot(aes(x = pre_mean, y = post_mean, color = imputation_reliable)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#4DAF4A", `FALSE` = "#E41A1C"),
                     labels = c(`TRUE` = "Reliable", `FALSE` = "Unreliable"),
                     name = NULL) +
  labs(title = "B. MNAR imputation audit",
       x = "Pre-imputation mean", y = "Post-imputation mean") +
  THM + theme(legend.position = "top")

# Panel C — Cohen's d histogram for MNAR proteins
d_median <- median(mnar_audit$effect_d, na.rm = TRUE)
d_iqr    <- IQR(mnar_audit$effect_d, na.rm = TRUE)

pC3 <- mnar_audit |>
  ggplot(aes(x = effect_d)) +
  geom_histogram(binwidth = 0.1, fill = "#984EA3", color = "white", alpha = 0.7) +
  geom_vline(xintercept = d_median, linetype = "dashed", color = "black") +
  annotate("text", x = d_median, y = Inf, vjust = 2, hjust = -0.1, size = 3.2,
           label = sprintf("median = %.2f\nIQR = %.2f", d_median, d_iqr)) +
  labs(title = "C. Imputation shift (MNAR proteins)",
       x = "Cohen's d (pre vs post)", y = "Count") +
  THM

imp_sub <- sprintf("missForest | %s values imputed | %d unreliable (>50%% missing)",
                   comma(n_total_miss), n_unreliable)

page2 <- pA3 / (pB3 | pC3) +
  plot_annotation(
    title = "Imputation Quality",
    subtitle = imp_sub,
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

pdf(file.path(REPORT_DIR, "02_imputation_report.pdf"), width = 14, height = 10)
print(page1)
print(page2)
dev.off()
cat("Wrote:", file.path(REPORT_DIR, "02_imputation_report.pdf"), "\n")
