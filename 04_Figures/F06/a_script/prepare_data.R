# F06 Within-Group Training-vs-Acute Concordance — Data Preparation
# Runs fGSEA for all state, within-group, and interaction contrasts.
# Caches results to c_data/ for downstream panels.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
})

set.seed(42)

DAT <- "04_Figures/F06/c_data"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# -- Load DEP results ---------------------------------------------------------

dep_df <- read_csv("03_DEP/a_non_imputed/c_data/03_combined_results.csv", show_col_types = FALSE)
message(sprintf("DEP: %d proteins x %d columns", nrow(dep_df), ncol(dep_df)))

# -- Build pathway collection --------------------------------------------------

pw <- build_pathway_collection(min_size = 15, max_size = 500)

# -- fGSEA for canonical figure contrasts -------------------------------------

contrasts <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_HR", "Training_LR", "Acute_HR", "Acute_LR",
  "Training_Interaction", "Acute_Interaction"
)
fgsea_list <- list()

for (ctr in contrasts) {
  t_col <- paste0("t_", ctr)
  if (!t_col %in% names(dep_df)) {
    message(sprintf("  Skipping %s — column %s not found", ctr, t_col))
    next
  }

  ranks <- setNames(dep_df[[t_col]], dep_df$gene)
  ranks <- ranks[!is.na(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]  # keep first occurrence
  ranks <- sort(ranks, decreasing = TRUE)

  res <- fgsea::fgseaMultilevel(
    pathways    = pw,
    stats       = ranks,
    minSize     = 15,
    maxSize     = 500,
    nPermSimple = 10000,
    eps         = 0
  )
  res <- as.data.frame(res)
  res$database <- classify_database(res$pathway)
  res$contrast <- ctr

  # Per-database BH correction
  res <- res %>%
    group_by(database) %>%
    mutate(padj = p.adjust(pval, method = "BH")) %>%
    ungroup()

  fgsea_list[[ctr]] <- tibble::as_tibble(res)
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  message(sprintf("  %s: %d pathways tested, %d sig (padj < 0.05)",
                  ctr, nrow(res), n_sig))
}

fgsea_all <- bind_rows(fgsea_list)

# Flatten leadingEdge to string for CSV export
fgsea_export <- fgsea_all %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  select(pathway, contrast, pval, padj, NES, ES, size, database, leadingEdge)

write_csv(fgsea_export, file.path(DAT, "fgsea_cache.csv"))

# GO:BP-only export
fgsea_gobp <- fgsea_export %>% filter(database == "GO:BP")
write_csv(fgsea_gobp, file.path(DAT, "fgsea_gobp.csv"))

message(sprintf("fGSEA cache: %d rows (%d GO:BP)", nrow(fgsea_export), nrow(fgsea_gobp)))

# -- Concordance stats for both pairs -----------------------------------------

pairs <- list(
  HR = list(x = "Training_HR", y = "Acute_HR",
            label_x = "Training", label_y = "Acute"),
  LR = list(x = "Training_LR", y = "Acute_LR",
            label_x = "Training", label_y = "Acute")
)

stats_list <- list()
for (pair_name in names(pairs)) {
  p <- pairs[[pair_name]]

  lfc_x <- dep_df[[paste0("logFC_", p$x)]]
  lfc_y <- dep_df[[paste0("logFC_", p$y)]]
  pi_x  <- dep_df[[paste0("pi_score_", p$x)]]
  pi_y  <- dep_df[[paste0("pi_score_", p$y)]]

  st <- concordance_stats(lfc_x, lfc_y, pi_x, pi_y, threshold = 0.05)

  stats_list[[pair_name]] <- tibble(
    pair            = pair_name,
    contrast_x      = p$x,
    contrast_y      = p$y,
    n_all           = st$n_all,
    pearson_all     = round(st$pearson_all, 4),
    pearson_ci_lo   = round(st$pearson_ci["lo"], 4),
    pearson_ci_hi   = round(st$pearson_ci["hi"], 4),
    spearman_all    = round(st$spearman_all, 4),
    spearman_ci_lo  = round(st$spearman_ci["lo"], 4),
    spearman_ci_hi  = round(st$spearman_ci["hi"], 4),
    sign_concordance_pct = round(st$sign_concordance_all, 1),
    n_sig           = if (!is.null(st$n_sig)) st$n_sig else NA_integer_,
    pearson_sig     = if (!is.null(st$pearson_sig)) round(st$pearson_sig, 4) else NA_real_,
    spearman_sig    = if (!is.null(st$spearman_sig)) round(st$spearman_sig, 4) else NA_real_,
    sign_conc_sig   = if (!is.null(st$sign_concordance_sig)) round(st$sign_concordance_sig, 1) else NA_real_
  )
}

conc_stats <- bind_rows(stats_list)
write_csv(conc_stats, file.path(DAT, "within_group_concordance_stats.csv"))

message("prepare_data.R done")
message(sprintf("  HR: r=%.3f, rho=%.3f, sign=%.1f%%",
                conc_stats$pearson_all[1], conc_stats$spearman_all[1],
                conc_stats$sign_concordance_pct[1]))
message(sprintf("  LR: r=%.3f, rho=%.3f, sign=%.1f%%",
                conc_stats$pearson_all[2], conc_stats$spearman_all[2],
                conc_stats$sign_concordance_pct[2]))
