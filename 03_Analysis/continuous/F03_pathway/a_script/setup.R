# F03 setup, continuous tree: read the DEP fit from upstream c_data and
# compute the fgsea enrichment for the two pooled contrasts (moderated-t
# ranks vs Hallmark, KEGG, Reactome, GO:BP/CC/MF, GO Slim), seeded. Mirrors
# categorical/F03_pathway/a_script/setup.R exactly except the contrast set
# -- this keeps the same full gene-set collection categorical's F03 uses
# (broader than 02_Pathways' Hallmark+GO-Slim-only set), matching that
# figure's own established scope rather than the feature layer's.
# Provides: dep, fg, pw, CONTRASTS, RPT_DIR, DAT_DIR. Writes nothing.

pacman::p_load(here, dplyr, tidyr, readr, tibble, fgsea, msigdbr, openxlsx)
source(here("functions", "shared_style.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("03_Features", "contrasts.R"))

RPT_DIR <- here("03_Analysis", "continuous", "F03_pathway", "b_reports")
DAT_DIR <- here("03_Analysis", "continuous", "F03_pathway", "c_data")

CONTRASTS <- trimws(sub("=.*$", "", POOLED_CONTRASTS))

dep <- read_csv(
  here(
    "03_Analysis", "continuous", "01_Proteins", "c_data",
    "02_combined_results.csv"
  ),
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
  res$leadingEdge <- vapply(
    res$leadingEdge, function(x) paste(x, collapse = ";"), character(1)
  )
  res
}) |>
  bind_rows()

if (!exists("F03_PANELS")) F03_PANELS <- list()
if (!exists("F03_REPORTS")) F03_REPORTS <- list()
if (!exists("F03_SUPP")) F03_SUPP <- list()
if (!exists("F03_TABLES")) F03_TABLES <- list()
