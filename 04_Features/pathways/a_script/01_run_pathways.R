# Pathway level of the feature layer: per-sample singscore ranks, then the same
# nine contrasts stage 03 fits on proteins. The scores are the feature space
# F05 and F06 predict from, so the cache written here is the one they read.
#
# Foroutan et al. 2018, BMC Bioinformatics 19:404 -- singscore
#
# Coverage travels with every row. A set can pass the 15-500 annotated-size
# filter with a tenth of its members measured, and a score built on 11 of 200
# proteins is not a readout of that pathway.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "feature_contrasts.R"))
source(here("functions", "shared_singscore.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "pred_features.R"))

out_dir <- here("04_Features", "pathways", "c_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(here("04_Features", "pathways", "b_reports"),
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

# The same sets tested a second way. singscore collapses a set to one score per
# sample and then tests that score; fry tests the member proteins directly,
# rotating residuals so inter-gene correlation is carried rather than assumed
# away. fgsea (F03) assumes it away, which is what makes the three-way
# comparison worth having. Complete data is required, so this runs on the
# missForest arm.
expr <- pred_gene_expression(readRDS(pred_paths()$dalist))
parts <- feature_design(expr)
fry_res <- pathway_fry(
  expr, gene_sets, parts$design, parts$contrasts,
  block = parts$block, correlation = parts$correlation
)

fry_summary <- fry_res |>
  group_by(.data$contrast) |>
  summarise(
    n_tested = dplyr::n(), n_nominal = sum(.data$p < 0.05),
    n_fdr = sum(.data$fdr < 0.05), min_fdr = min(.data$fdr), .groups = "drop"
  )
cat("\nfry (rotation test, subject-blocked) over the same sets:\n")
print(as.data.frame(fry_summary))

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
    fry = fry_res, fry_summary = fry_summary
  ),
  file.path(out_dir, "pathway_contrasts.xlsx")
)
