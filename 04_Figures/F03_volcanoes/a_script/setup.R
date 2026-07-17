# F03 setup: read the DEP fit from upstream c_data and compute the fgsea enrichment
# for all 9 contrasts (moderated-t ranks vs Hallmark, KEGG, Reactome, GO:BP/CC/MF, GO
# Slim), seeded. The cache is RAW - no redundancy collapse - because the EnrichmentMap
# dedup is a display decision and belongs at ring-draw time, not in the cache.
# Provides: dep, fg, pw, CONTRASTS, RPT_DIR, DAT_DIR. Writes nothing.

pacman::p_load(here, dplyr, tidyr, readr, tibble, fgsea, msigdbr, openxlsx)
source(here("functions", "shared_style.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("03_DEP", "contrasts.R"))

RPT_DIR <- here("04_Figures", "F03_volcanoes", "b_reports")
DAT_DIR <- here("04_Figures", "F03_volcanoes", "c_data")

CONTRASTS <- sub(" =.*$", "", HRVLR_CONTRASTS)

dep <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)

pw <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
)
set.seed(42)
fg <- lapply(CONTRASTS, function(ct) {
  d <- tibble(gene = dep$gene, t = dep[[paste0("t_", ct)]]) |>
    filter(!is.na(gene), !is.na(t)) |>
    distinct(gene, .keep_all = TRUE)
  ranks <- sort(setNames(d$t, d$gene), decreasing = TRUE)
  res <- run_fgsea(ranks, pw)
  res$contrast <- ct
  res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  res
}) |>
  bind_rows()


if (!exists("F03_PANELS")) F03_PANELS <- list()
if (!exists("F03_REPORTS")) F03_REPORTS <- list()
if (!exists("F03_SUPP")) F03_SUPP <- list()
if (!exists("F03_TABLES")) F03_TABLES <- list()
