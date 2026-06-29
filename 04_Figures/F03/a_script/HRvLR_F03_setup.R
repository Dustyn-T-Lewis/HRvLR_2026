# HRvLR_F03_setup.R - Figure 3 (enrichVolcano ring-volcanoes).
# Computes the fgsea enrichment cache ONCE for all 9 DEP contrasts (moderated-t
# ranks vs Hallmark, KEGG, Reactome, GO:BP, and GO Slim). It feeds the volcanoes
# and the F04/F05 NES scatters - computed once, read by all three.
# Provides: dep (combined DEP results), fg (fgsea cache), CONTRASTS, RPT_DIR,
# DAT_DIR + style.R / pathway_utils.R exports. The cache lives in the one
# workbook F03_source_data.xlsx; delete it to recompute.

pacman::p_load(here, dplyr, tidyr, readr, tibble, fgsea, msigdbr, openxlsx)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))
set.seed(42)

RPT_DIR <- here("04_Figures", "F03", "b_reports")
DAT_DIR <- here("04_Figures", "F03", "c_data")
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c(
  "Training_HR", "Training_LR", "Acute_HR", "Acute_LR",
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)

dep <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)

cache_xlsx <- file.path(DAT_DIR, "F03_source_data.xlsx")
if (file.exists(cache_xlsx)) {
  fg <- as_tibble(read.xlsx(cache_xlsx, sheet = "fgsea_all"))
} else {
  pw <- build_pathway_collection(
    min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
  )
  fg <- lapply(CONTRASTS, function(ct) {
    d <- tibble(gene = dep$gene, t = dep[[paste0("t_", ct)]]) |>
      filter(!is.na(gene), !is.na(t)) |>
      distinct(gene, .keep_all = TRUE)
    ranks <- sort(setNames(d$t, d$gene), decreasing = TRUE)
    res <- run_fgsea_deduplicated(ranks, pw)
    res$contrast <- ct
    res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    res
  }) |>
    bind_rows()
  sig <- fg |>
    filter(!is.na(padj), padj < 0.05, is.finite(NES)) |>
    transmute(
      contrast, database,
      pathway = clean_pathway_name(pathway),
      NES, padj, size, leading_edge = leadingEdge
    ) |>
    arrange(contrast, padj)
  wb <- createWorkbook()
  addWorksheet(wb, "fgsea_all")
  writeData(wb, "fgsea_all", fg)
  addWorksheet(wb, "fgsea_significant")
  writeData(wb, "fgsea_significant", sig)
  saveWorkbook(wb, cache_xlsx, overwrite = TRUE)
}

cat(sprintf("F03 setup: %d contrasts, fgsea cache %d rows\n", length(CONTRASTS), nrow(fg)))
