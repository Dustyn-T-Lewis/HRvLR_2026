suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_assoc_leaf.R"))
}))

bundle <- pred_load()
grid <- expand.grid(
  level = SWEEP_LEVELS, config = SWEEP_CONFIGS, method = ASSOC_METHODS,
  stringsAsFactors = FALSE
)

tidy_all <- dplyr::bind_rows(Map(function(l, cfg, m) {
  message(sprintf("[assoc] %s · %s · %s", l, cfg, m))
  run_assoc_leaf(bundle, l, cfg, m)
}, grid$level, grid$config, grid$method))

summary_all <- assoc_cell_summary(tidy_all)

# Nominal, not BH: F04-F06 report a raw p against a screen-size denominator,
# and a q here would imply an independence the cells do not have. The `bh`
# column stays in the leaf sheets; no panel counts it.
root <- sweep_root_dir("F04_association")
dir.create(file.path(root, "c_data"), recursive = TRUE, showWarnings = FALSE)
write_sweep_workbook(
  file.path(root, "c_data", "results.xlsx"),
  list(
    cell_summary = summary_all,
    hits = dplyr::filter(tidy_all, .data$p < 0.05)
  )
)
message(sprintf(
  "association done: %d cells, %d nominal hits",
  nrow(summary_all), sum(summary_all$n_nominal)
))
