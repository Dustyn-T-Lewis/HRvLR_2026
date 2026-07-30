# Module level of the feature layer: the signed WGCNA eigengenes from
# 01_run_modules.R, put through the same nine contrasts stage 03 fits on
# proteins. Nothing upstream tests eigengenes against the design, so this is the
# only place the module answer to the responder question is computed.
#
# Langfelder & Horvath 2008, BMC Bioinformatics 9:559 -- WGCNA
#
# Two limits to state wherever these numbers appear. Module detection is
# unsupervised, so it never saw the HR/LR labels, but it did see all 45 samples,
# and eBayes has 11 features to borrow variance across, which is close to no
# moderation at all.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "feature_contrasts.R"))

out_dir <- here("04_Features", "modules", "c_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eigengenes <- module_matrix()
results <- fit_feature_contrasts(eigengenes)

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
