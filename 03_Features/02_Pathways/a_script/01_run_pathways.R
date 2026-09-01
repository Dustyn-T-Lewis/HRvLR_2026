# Pathway level of the feature layer: per-sample singscore ranks, tested
# through the same estimator the protein level fits, plus a second,
# independent per-contrast test (fgsea) over the same sets.
#
# Both contrast families run here off one set of scores. The nine
# HR-vs-LR contrasts and the two pooled Training/Acute contrasts differ
# only in their design matrix, so fitting them separately and binding on a
# family column keeps one set of sets, one cache and one workbook where
# there used to be two directories that could drift apart.
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

out_dir <- here("03_Features", "02_Pathways", "c_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(here("03_Features", "02_Pathways", "b_reports"),
  recursive = TRUE, showWarnings = FALSE
)

# Fit one feature matrix through both designs and tag which is which.
fit_both <- function(mat) {
  bind_rows(
    fit_feature_contrasts(mat) |> mutate(family = "categorical"),
    fit_feature_contrasts_pooled(mat) |> mutate(family = "pooled")
  )
}

collection <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
)
gene_sets <- collection[
  classify_database(names(collection)) %in% c("Hallmark", "GO Slim")
]

scores <- pathway_matrix()
coverage <- pathway_coverage(gene_sets, detected_genes()) |>
  filter(.data$feature %in% rownames(scores))

results <- fit_both(scores) |>
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
expr_fit <- fit_both(expr)

cells <- distinct(expr_fit, .data$family, .data$contrast)
fgsea_res <- bind_rows(lapply(seq_len(nrow(cells)), function(i) {
  rows <- filter(
    expr_fit,
    .data$family == cells$family[i], .data$contrast == cells$contrast[i]
  )
  run_fgsea(sort(setNames(rows$t, rows$feature)), gene_sets) |>
    mutate(family = cells$family[i], contrast = cells$contrast[i])
}))

fgsea_summary <- fgsea_res |>
  group_by(.data$family, .data$contrast) |>
  summarise(
    n_tested = dplyr::n(), n_nominal = sum(.data$pval < 0.05),
    n_fdr = sum(.data$padj < 0.05), min_fdr = min(.data$padj), .groups = "drop"
  )
cat("\nfgsea (preranked) over the same sets:\n")
print(as.data.frame(fgsea_summary))

summary_tbl <- results |>
  group_by(.data$family, .data$contrast) |>
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
      .data$family, .data$contrast, .data$feature, .data$n_detected,
      .data$n_annotated,
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
