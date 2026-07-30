# Protein level of the feature layer. Stage 03 already fits the nine contrasts
# on the non-imputed cycloess matrix through proteoDA, so this reshapes that fit
# into the long form the other two levels write, and re-checks it against a
# refit.
#
# Ritchie et al. 2015, Nucleic Acids Res 43(7):e47 -- limma
# Smyth, Michaud & Scott 2005, Bioinformatics 21(9):2067 -- duplicateCorrelation

pacman::p_load(here, dplyr, purrr, openxlsx)
source(here("functions", "feature_contrasts.R"))

out_dir <- here("04_Features", "proteins", "c_data")
report_dir <- here("04_Features", "proteins", "b_reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

book <- here("03_DEP", "a_non_imputed", "c_data", "05_results.xlsx")
contrast_names <- trimws(sub("=.*$", "", HRVLR_CONTRASTS))

results <- map_dfr(contrast_names, function(ct) {
  read.xlsx(book, ct) |>
    transmute(
      contrast = ct, feature = .data$uniprot_id, gene = .data$gene,
      logFC = .data$logFC, p = .data$P.Value, bh = .data$adj.P.Val
    )
})

equivalence <- verify_protein_equivalence()
if (!all(equivalence$equivalent)) {
  stop("stage 03 workbook and a refit through feature_contrasts disagree")
}

summary_tbl <- results |>
  group_by(.data$contrast) |>
  summarise(
    n_tested = sum(!is.na(.data$p)),
    n_nominal = sum(.data$p < 0.05, na.rm = TRUE),
    n_bh = sum(.data$bh < 0.05, na.rm = TRUE),
    min_bh = min(.data$bh, na.rm = TRUE), .groups = "drop"
  )
print(as.data.frame(summary_tbl))

write.xlsx(
  list(contrasts = results, summary = summary_tbl, equivalence = equivalence),
  file.path(out_dir, "protein_contrasts.xlsx")
)
