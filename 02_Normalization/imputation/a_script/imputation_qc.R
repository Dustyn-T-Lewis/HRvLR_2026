#!/usr/bin/env Rscript
# Imputation QC: the four arms compared against the reported non-imputed fit.
#
# This lives with the imputation it documents, but it reads 03_DEP/b_imputed, because the only
# honest way to compare imputers is by what they do to the test downstream. Run it after stage 03.
#
# On BH the null holds under every imputer that does not impute inside the tested factor:
# missForest (MAR), MsCoreUtils (hybrid) and Perseus (MNAR) each return zero BH<0.10 hits in
# all five HR-vs-LR and interaction contrasts. imp4p returns 117, because impute.mle fits a
# separate EM within each Group_Time cell and the contrasts then test among those same cells;
# re-imputing within random cells of equal size collapses it back to zero (02_imp4p_circularity.R).
#
# Do NOT rank robustness on pi-counts. pi = p^|log2FC| exponentiates the fold change, which is
# exactly what imputation distorts: Perseus has the MOST pi-hits of any arm and the SAME BH count
# as the non-imputed fit. missForest stays canonical downstream.

pacman::p_load(here, dplyr, tidyr, readr, purrr, ggplot2, patchwork, openxlsx)
source(here("04_Figures", "functions", "shared_style.R"))

RPT_DIR <- here("02_Normalization", "imputation", "b_reports")
DAT_DIR <- here("02_Normalization", "imputation", "c_data")

# All nine, ordered so the four within-arm contrasts come first and the five the paper reports as
# null are grouped on the right. Showing only a subset was how the old figure managed to omit four
# of the five contrasts its claim depended on.
within_arm <- c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR")
null_contrasts <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)
key_contrasts <- c(within_arm, null_contrasts)
strata_levels <- c("0 (imputation is a no-op)", "1-9", ">=10")
methods <- c("missforest", "imp4p", "mscoreutils", "perseus")
method_lab <- c(
  missforest = "missForest", imp4p = "imp4p",
  mscoreutils = "MsCoreUtils", perseus = "Perseus"
)

method_levels <- c("Non-imputed", unname(method_lab))
relabel <- function(d) {
  d |>
    filter(contrast %in% key_contrasts) |>
    mutate(
      method = factor(c(non_imputed = "Non-imputed", method_lab)[method], levels = method_levels),
      contrast = factor(contrast, levels = key_contrasts)
    )
}

# BH beside Pi, from the table stage 03 writes. BH is the metric the headline rests on; Pi is
# shown next to it because it is what the earlier version of this figure reported on its own, and
# on its own it misleads. Pi = p^|log2FC| exponentiates the fold change, which is exactly what
# imputation distorts: Perseus carries the most Pi-hits of any arm and the same BH count as the
# non-imputed fit.
sens <- relabel(read_csv(
  here("03_DEP", "b_imputed", "c_data", "sensitivity_bh_vs_pi.csv"),
  show_col_types = FALSE
))

p_bh <- ggplot(sens, aes(contrast, n_BH_10, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(
    title = "BH-significant proteins (q < 0.10): the metric the result rests on",
    subtitle = "Zero under every imputer except imp4p, which imputes inside the Group_Time cells the contrasts test",
    x = NULL, y = "BH q < 0.10", tag = "A"
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")

p_pi <- ggplot(sens, aes(contrast, n_pi, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(
    title = "Pi-selected proteins: not comparable across arms",
    subtitle = "Pi exponentiates |log2FC|, the quantity imputation distorts. Perseus tops this panel and is null in panel A.",
    x = NULL, y = "Pi < 0.05", tag = "B"
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

# Concordance, split by how much of a protein was actually missing. Unstratified this metric is
# worthless: the ~944 fully observed proteins are a no-op for every imputer and force rho to 1,
# which is why it once ranked imp4p second best while imp4p was destroying the null.
conc <- read_csv(
  here("03_DEP", "b_imputed", "c_data", "logfc_concordance.csv"),
  show_col_types = FALSE
) |>
  relabel() |>
  mutate(stratum = factor(stratum, levels = strata_levels)) |>
  group_by(method, stratum) |>
  summarise(rho = median(rho, na.rm = TRUE), .groups = "drop")

p_conc <- ggplot(conc, aes(stratum, method, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = FIG_GEOM_TEXT, fontface = "bold") +
  scale_fill_gradient(low = "#FEE0D2", high = "#67000D", name = "median rho\nvs non-imputed") +
  labs(
    title = "logFC concordance, by how many values were missing",
    subtitle = "A rank correlation cannot see imp4p's failure, which is a variance collapse. Diagnostic, not a safeguard.",
    x = "missing values per protein (of 45)", y = NULL, tag = "C"
  ) +
  FIG_THEME

composite <- (p_bh / p_pi / p_conc) +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(
    title = "Imputation QC. The null holds under every imputer that does not impute inside the tested factor",
    subtitle = "Non-imputed reported arm vs four imputation methods. missForest is the canonical downstream matrix.",
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
    )
  )

ggsave(file.path(RPT_DIR, "imputation_qc.png"), composite,
  width = 220, height = 240, units = "mm", dpi = 200, bg = "white"
)
ggsave(file.path(RPT_DIR, "imputation_qc.pdf"), composite,
  width = 220, height = 240, units = "mm", device = PDF_DEVICE, bg = "white"
)

overview <- data.frame(
  sheet = c("bh_vs_pi", "logfc_concordance"),
  description = c(
    "Panels A-B: BH q<0.10 and Pi<0.05 protein counts per method per key contrast",
    "Panel C: median Spearman rho of each imputed arm's per-protein logFC vs the non-imputed reference, split by per-protein missingness"
  ),
  stringsAsFactors = FALSE
)
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
addWorksheet(wb, "bh_vs_pi")
writeData(wb, "bh_vs_pi", as.data.frame(sens))
addWorksheet(wb, "logfc_concordance")
writeData(wb, "logfc_concordance", as.data.frame(conc))
saveWorkbook(wb, file.path(DAT_DIR, "imputation_qc_source_data.xlsx"), overwrite = TRUE)
cat("imputation QC composite + workbook saved to", RPT_DIR, "\n")
