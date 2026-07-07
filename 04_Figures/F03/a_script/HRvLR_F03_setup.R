# HRvLR_F03_setup.R - Figure 3 (enrichVolcano ring-volcanoes).
# Computes the fgsea enrichment for all 9 DEP contrasts (moderated-t ranks vs
# Hallmark, KEGG, Reactome, GO:BP, and GO Slim) fresh every run, seeded for
# reproducibility. It feeds the volcanoes and the F04/F05 NES scatters.
# Provides: dep (combined DEP results), fg (fgsea cache), CONTRASTS, RPT_DIR,
# DAT_DIR + style.R / pathway_utils.R exports. The workbook
# F03_source_data.xlsx is overwritten fresh every run.

pacman::p_load(here, dplyr, tidyr, readr, tibble, fgsea, msigdbr, openxlsx)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))

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

pw <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
)
set.seed(42)
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

overview <- data.frame(
  sheet = c("fgsea_all", "fgsea_significant"),
  description = c(
    "All fgsea enrichment results across the 9 DEP contrasts; source for the ring-volcano panels",
    "Significant subset (padj < 0.05), cleaned pathway names, arranged by contrast and padj"
  ),
  stringsAsFactors = FALSE
)
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
addWorksheet(wb, "fgsea_all")
writeData(wb, "fgsea_all", fg)
addWorksheet(wb, "fgsea_significant")
writeData(wb, "fgsea_significant", sig)
saveWorkbook(wb, cache_xlsx, overwrite = TRUE)

cat(sprintf("F03 setup: %d contrasts, fgsea cache %d rows\n", length(CONTRASTS), nrow(fg)))
