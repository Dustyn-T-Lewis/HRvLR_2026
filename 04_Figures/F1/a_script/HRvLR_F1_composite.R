# Figure 1 — Composite Assembly + Supplementary Excel Workbook ----
#
# Sources all per-panel scripts, assembles with patchwork into a single
# Figure_1.pdf/png, and exports F1_supplementary.xlsx.
#
# Panels:
#   A — CV% violins (inter-individual variability)
#   B — |logFC| histograms (effect size distributions)
#   C — PCA biplot + PERMANOVA + legend key
#   D — DEPs per contrast (pseudo-log stacked bar)
#   E — UpSet plot (dual bar + dot matrix + legend key)
#   F — fGSEA grouped bar chart (pathway enrichment)
#   S1.7 — Pi-score distributions (supplementary)

# === 0. Source all panels (setup loaded by panel_A) ==========================

source("04_Figures/F1/a_script/panel_A.R")     # -> pA
source("04_Figures/F1/a_script/panel_B.R")     # -> pB
source("04_Figures/F1/a_script/panel_C.R")     # -> pC_group, pC_time, pC_combined
source("04_Figures/F1/a_script/panel_D.R")     # -> pD (+ sig_sets, dir_map for E)
source("04_Figures/F1/a_script/panel_E.R")     # -> pE_bars, pE_dots, pKeys, pE_combined
source("04_Figures/F1/a_script/panel_F.R")     # -> pF
source("04_Figures/F1/a_script/supp_S1_7.R")   # -> s17

message("All F1 panel scripts sourced — assembling composite...")

# === 1. Composite figure (patchwork) =========================================

left_col <- (pA / pB / pC_combined) +
  plot_layout(heights = c(0.26, 0.34, 0.40))

mid_col <- (pD / pE_combined) +
  plot_layout(heights = c(0.28, 0.72))

fig1 <- (left_col | mid_col | pF) +
  plot_layout(widths = c(0.37, 0.38, 0.25)) +
  plot_annotation(
    title = "Proteomic Landscape of Resistance Training Response in High and Low Responders",
    subtitle = "Data quality, effect magnitude, contrast overlap, and functional enrichment across the 2 \u00d7 3 factorial design.",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT_DIR, "Figure_1.pdf"), fig1,
       width = 400, height = 260, units = "mm", device = PDF_DEVICE)
ggsave(file.path(RPT_DIR, "Figure_1.png"), fig1,
       width = 400, height = 260, units = "mm", dpi = 300)
cat("Saved Figure_1.pdf and Figure_1.png\n")

# === 2. Supplementary Excel workbook =========================================

library(openxlsx)

# Sheet A: variability
sheet_A <- cv_df |>
  left_join(cv_med, by = "group") |>
  rename(cv_pct = cv, group_median_cv = med)

# Sheet B: effect sizes
sheet_B <- lfc_stats |>
  select(contrast, med_lfc, med_abs_lfc, ci_lo, ci_hi, n_above_05)

# Sheet C: PCA scores
sheet_C <- pca_df |>
  select(sample_id, PC1, PC2, Group, Timepoint, group)

# Sheet D: DEP fractions
sheet_D <- frac_df |>
  select(contrast, threshold, n, pct)

# Sheet E: upset intersections
upset_tbl <- tibble(
  intersection_id = seq_len(n_int),
  comb_code       = comb_names_ord,
  up_count        = up_ord,
  down_count      = down_ord,
  total           = up_ord + down_ord
)
for (j in seq_along(set_order_ch)) {
  upset_tbl[[set_order_ch[j]]] <- vapply(comb_names_ord, function(code) {
    as.logical(as.integer(strsplit(code, "")[[1]])[j])
  }, logical(1))
}
sheet_E <- upset_tbl

# Sheet F: fGSEA results
sheet_F <- fgsea_export

# Sheet meta
sheet_meta <- tibble(
  item = c("figure_title", "description", "contrasts", "statistical_methods",
           "sample_sizes", "significance_threshold", "software", "date_generated"),
  value = c(
    "Figure 1: Proteomic Landscape of Resistance Training Response in High and Low Responders",
    "Data quality (CV%), effect magnitude (|logFC|), PCA, DEP counts, contrast overlap (UpSet), and pathway enrichment (fGSEA).",
    paste(CONTRASTS, collapse = ", "),
    "Wilcoxon rank-sum (CV%), KS test (logFC), PERMANOVA (PCA), pi-score (DEPs), fGSEA multilevel (pathways), Jaccard dedup (0.5)",
    sprintf("%d samples, %d proteins (filtered), %d DEP results",
            nrow(meta), nrow(norm_df), nrow(dep_df)),
    "Pi-score < 0.05 (DEPs), padj < 0.05 (fGSEA)",
    "R, limma, fgsea, patchwork, openxlsx",
    format(Sys.Date(), "%Y-%m-%d")
  )
)

wb <- createWorkbook()
sheets <- list(A_variability = sheet_A, B_effect_sizes = sheet_B,
               C_pca_scores = sheet_C, D_dep_fractions = sheet_D,
               E_upset_intersections = sheet_E, F_fgsea_results = sheet_F,
               `_metadata` = sheet_meta)
hs <- createStyle(textDecoration = "bold")
for (sn in names(sheets)) {
  addWorksheet(wb, sn)
  writeData(wb, sn, sheets[[sn]])
  addStyle(wb, sn, hs, rows = 1, cols = seq_len(20), gridExpand = TRUE)
  setColWidths(wb, sn, cols = seq_len(20), widths = "auto")
}

xlsx_path <- file.path(DAT_DIR, "F1_supplementary.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat(sprintf("Saved %s (%d sheets)\n", xlsx_path, length(sheets)))

cat("Figure 1 composite assembly complete.\n")
