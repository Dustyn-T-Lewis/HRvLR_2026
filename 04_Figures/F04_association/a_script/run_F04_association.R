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
  message(sprintf("[assoc] %s | %s | %s", l, cfg, m))
  run_assoc_leaf(bundle, l, cfg, m)
}, grid$level, grid$config, grid$method))

summary_all <- assoc_cell_summary(tidy_all)

root <- sweep_root_dir("F04_association")
dir.create(file.path(root, "c_data"), recursive = TRUE, showWarnings = FALSE)
write_sweep_workbook(
  file.path(root, "c_data", "results.xlsx"),
  list(
    cell_summary = summary_all,
    bh_hits = dplyr::filter(tidy_all, .data$bh < 0.05)
  )
)
message(sprintf(
  "association done: %d cells, %d BH survivors",
  nrow(summary_all), sum(summary_all$n_bh)
))
