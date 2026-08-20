# Pathway level of the feature layer: per-sample singscore ranks, tested
# through the same estimator the nine contrasts fit on proteins, plus a
# second, independent per-contrast test (fgsea) over the same sets. The
# singscore matrix cached here is the sample-level pathway feature
# F04_classification reads.
#
# Foroutan et al. 2018, BMC Bioinformatics 19:404 -- singscore
#
# Coverage travels with every row. A set can pass the 15-500 annotated-size
# filter with a tenth of its members measured, and a score built on 11 of 200
# proteins is not a readout of that pathway.
#
# fgsea has no rotation-test cross-check here (fry was dropped from every
# pathway/module call in this tree). fgsea permutes gene labels and assumes
# they vary independently, which they do not on shared biopsies -- this
# repo caught a concrete false positive this way before (an OxPhos call at
# padj 4e-17 that a rotation test rejected). Without that check, read a
# large fgsea hit count as what preranked GSEA alone can find, not as
# confirmed enrichment.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "feature_contrasts.R"))
source(here("functions", "shared_singscore.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "pred_features.R"))

out_dir <- here("03_Analysis", "categorical", "02_Pathways", "c_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(here("03_Analysis", "categorical", "02_Pathways", "b_reports"),
  recursive = TRUE, showWarnings = FALSE
)

collection <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
)
gene_sets <- collection[
  classify_database(names(collection)) %in% c("Hallmark", "GO Slim")
]

scores <- pathway_matrix()
coverage <- pathway_coverage(gene_sets, detected_genes()) |>
  filter(.data$feature %in% rownames(scores))

results <- fit_feature_contrasts(scores) |>
  left_join(coverage, by = "feature") |>
  mutate(
    database = classify_database(.data$feature),
    pi = pi_score(.data$p, .data$logFC)
  )

# The same sets tested a second way, per contrast: singscore collapses a set
# to one score per sample and tests that score through the shared estimator
# above; fgsea asks whether the set's own genes cluster at one end of that
# contrast's full ranked list (canonical fgseaMultilevel, run_fgsea() in
# shared_pathway_utils.R). Complete data is required, so this runs on the
# missForest arm, matching pred_gene_expression()'s gene-symbol collapse --
# the same matrix msigdbr's gene sets are keyed on.
expr <- pred_gene_expression(readRDS(pred_paths()$dalist))
expr_fit <- fit_feature_contrasts(expr)

fgsea_res <- bind_rows(lapply(unique(expr_fit$contrast), function(ct) {
  rows <- filter(expr_fit, .data$contrast == ct)
  run_fgsea(sort(setNames(rows$t, rows$feature)), gene_sets) |>
    mutate(contrast = ct)
}))

fgsea_summary <- fgsea_res |>
  group_by(.data$contrast) |>
  summarise(
    n_tested = dplyr::n(), n_nominal = sum(.data$pval < 0.05),
    n_fdr = sum(.data$padj < 0.05), min_fdr = min(.data$padj), .groups = "drop"
  )
cat("\nfgsea (preranked) over the same sets:\n")
print(as.data.frame(fgsea_summary))

summary_tbl <- results |>
  group_by(.data$contrast) |>
  summarise(
    n_tested = dplyr::n(), n_nominal = sum(.data$p < 0.05),
    n_bh = sum(.data$bh < 0.05), min_bh = min(.data$bh), .groups = "drop"
  )
print(as.data.frame(summary_tbl))

cat(sprintf(
  "\n%d of %d sets have under 25%% of their annotated members detected\n",
  sum(coverage$coverage < 0.25), nrow(coverage)
))
print(as.data.frame(
  results |>
    filter(.data$bh < 0.05) |>
    arrange(.data$bh) |>
    transmute(
      .data$contrast, .data$feature, .data$n_detected, .data$n_annotated,
      logFC = round(.data$logFC, 3), bh = signif(.data$bh, 3)
    )
))

write.xlsx(
  list(
    contrasts = results, summary = summary_tbl, coverage = coverage,
    fgsea = fgsea_res, fgsea_summary = fgsea_summary
  ),
  file.path(out_dir, "pathway_contrasts.xlsx")
)
