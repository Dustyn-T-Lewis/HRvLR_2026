# Module level of the feature layer, continuous tree: the signed WGCNA
# eigengenes from 01_run_modules.R, put through the same two pooled
# contrasts stage 03 fits on proteins. Mirrors
# categorical/03_WGCNA/a_script/02_run_module_contrasts.R -- nothing
# upstream tests eigengenes against the design, so this is the only place
# the module answer to the training/acute question is computed.
#
# Langfelder & Horvath 2008, BMC Bioinformatics 9:559 -- WGCNA
#
# Two limits to state wherever these numbers appear. Module detection is
# unsupervised and never saw Timepoint either, but it did see all 45
# samples, and eBayes has as many features as there are modules to borrow
# variance across, which is close to no moderation at all.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "feature_contrasts.R"))

out_dir <- here("03_Analysis", "continuous", "03_WGCNA", "c_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eigengenes <- module_matrix(tree = "continuous")
results <- fit_feature_contrasts_pooled(eigengenes)

summary_tbl <- results |>
  group_by(.data$contrast) |>
  summarise(
    n_tested = dplyr::n(), n_nominal = sum(.data$p < 0.05),
    n_bh = sum(.data$bh < 0.05), min_bh = min(.data$bh), .groups = "drop"
  )

cat(sprintf(
  "%d modules x %d samples | within-subject correlation %.3f\n",
  nrow(eigengenes), ncol(eigengenes), attr(results, "within_cor")
))
print(as.data.frame(summary_tbl))

write.xlsx(
  list(contrasts = results, summary = summary_tbl),
  file.path(out_dir, "module_contrasts.xlsx")
)
