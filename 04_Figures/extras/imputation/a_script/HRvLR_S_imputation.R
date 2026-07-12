#!/usr/bin/env Rscript
# Supplement S: imputation-method comparison on the key contrasts.
# Non-imputed reported arm vs missForest / imp4p / MsCoreUtils / Perseus. Shows the
# key effects hold across imputation choice. missForest stays canonical downstream.

pacman::p_load(here, dplyr, tidyr, readr, purrr, ggplot2, patchwork, openxlsx)
source(here("04_Figures", "functions", "style.R"))

RPT_DIR <- here("04_Figures", "extras", "imputation", "b_reports")
DAT_DIR <- here("04_Figures", "extras", "imputation", "c_data")

key_contrasts <- c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR", "Baseline_HRvLR")
methods <- c("missforest", "imp4p", "mscoreutils", "perseus")
method_lab <- c(
  missforest = "missForest", imp4p = "imp4p",
  mscoreutils = "MsCoreUtils", perseus = "Perseus"
)

concordance <- read_csv(
  here("03_DEP", "b_imputed", "c_data", "logfc_concordance.csv"),
  show_col_types = FALSE
) |>
  filter(contrast %in% key_contrasts, method %in% methods) |>
  mutate(
    method = factor(method_lab[method], levels = unname(method_lab)),
    contrast = factor(contrast, levels = key_contrasts)
  )

count_imp <- map_dfr(methods, function(m) {
  read_csv(
    here("03_DEP", "b_imputed", "c_data", m, "combined_results_pi.csv"),
    show_col_types = FALSE
  ) |>
    filter(contrast %in% key_contrasts) |>
    group_by(contrast) |>
    summarise(n_pi = sum(sig_pi != 0, na.rm = TRUE), .groups = "drop") |>
    mutate(method = m)
})
ni <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)
count_ni <- tibble(
  contrast = key_contrasts,
  n_pi = vapply(key_contrasts, function(ct) {
    sum(ni[[paste0("sig_pi_", ct)]] != 0, na.rm = TRUE)
  }, integer(1)),
  method = "non_imputed"
)
counts <- bind_rows(count_ni, count_imp) |>
  filter(contrast %in% key_contrasts) |>
  mutate(
    method = factor(
      c(non_imputed = "Non-imputed", method_lab)[method],
      levels = c("Non-imputed", unname(method_lab))
    ),
    contrast = factor(contrast, levels = key_contrasts)
  )

p_conc <- ggplot(concordance, aes(contrast, method, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = FIG_GEOM_TEXT, fontface = "bold") +
  scale_fill_gradient(
    low = "#FEE0D2", high = "#67000D", limits = c(min(0.9, min(concordance$rho)), 1),
    name = "Spearman ρ\nvs non-imputed"
  ) +
  labs(title = "logFC concordance with the non-imputed arm", x = NULL, y = NULL, tag = "A") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_counts <- ggplot(counts, aes(contrast, n_pi, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(title = "π-significant proteins per method", x = NULL, y = "π-significant", tag = "B") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

composite <- (p_conc / p_counts) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Supplement S. Imputation-method comparison on the key contrasts",
    subtitle = "Non-imputed reported arm vs four imputation methods. missForest is the canonical downstream matrix.",
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
    )
  )

ggsave(file.path(RPT_DIR, "S_imputation_composite.png"), composite,
  width = 220, height = 240, units = "mm", dpi = 200, bg = "white"
)
ggsave(file.path(RPT_DIR, "S_imputation_composite.pdf"), composite,
  width = 220, height = 240, units = "mm", device = PDF_DEVICE, bg = "white"
)

overview <- data.frame(
  sheet = c("logfc_concordance", "pi_counts"),
  description = c(
    "Panel A: Spearman rho of each imputed arm's per-protein logFC vs the non-imputed reference, key contrasts",
    "Panel B: pi-significant protein count per method per key contrast"
  ),
  stringsAsFactors = FALSE
)
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
addWorksheet(wb, "logfc_concordance")
writeData(wb, "logfc_concordance", as.data.frame(concordance))
addWorksheet(wb, "pi_counts")
writeData(wb, "pi_counts", as.data.frame(counts))
saveWorkbook(wb, file.path(DAT_DIR, "S_imputation_source_data.xlsx"), overwrite = TRUE)
cat("S_imputation composite + workbook saved to", RPT_DIR, "\n")
