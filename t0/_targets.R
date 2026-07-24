# Single entry point. Run from the repo root so renv activates:
#   targets::tar_make(script = "t0/_targets.R", store = "t0/_targets")
# Set R_CONFIG_ACTIVE=quick for a fast validation slice. Every step is a pure function from t0/R/.

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "dplyr", "tidyr", "tibble", "purrr", "stringr", "readr", "here", "config",
    "rlang", "recipes", "parsnip", "workflows", "workflowsets", "tune", "rsample",
    "yardstick", "matrixStats", "singscore", "WGCNA", "msigdbr", "janitor",
    "ggplot2", "withr", "sessioninfo", "methods", "stats"
  ),
  seed = 20260723L
)

tar_source(here::here("t0", "R"))

write_out <- function(x, sub, name, cfg) {
  dir <- here::here(cfg$paths$outputs, sub)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, name)
  readr::write_csv(x, path)
  path
}

write_session_info <- function(cfg) {
  dir <- here::here(cfg$paths$outputs)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, "session_info.txt")
  writeLines(utils::capture.output(sessioninfo::session_info()), path)
  path
}

# The Quarto report is optional: without the quarto R package and CLI the rest of the DAG still
# builds. Install both to render it.
quarto_ready <- requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))

targets <- list(
  tar_target(cfg, config::get(file = here::here("t0", "config.yml"))),
  tar_target(proteins, load_proteins(cfg)),
  tar_target(meta, load_meta(cfg)),
  tar_target(valid, validate_inputs(proteins, meta, cfg)),
  tar_target(gene_sets, load_gene_sets(cfg)),
  tar_target(gmap, gene_map(proteins, cfg)),
  tar_target(qc_miss, qc_missingness(proteins, meta)),
  tar_target(qc_dist, qc_distributions(proteins, meta)),
  tar_target(qc_pair, qc_pairing(meta, cfg)),
  tar_target(qc_out, qc_outliers(proteins, cfg)),
  tar_target(qc_out_csv, write_out(qc_out, "qc", "outlier_flags.csv", cfg), format = "file"),
  tar_target(qc_pair_csv, write_out(qc_pair, "qc", "pairing.csv", cfg), format = "file"),
  tar_target(unsupervised, run_unsupervised(proteins, meta, gene_sets, gmap, cfg)),
  tar_target(unsup_csv, write_out(unsupervised, "results", "unsupervised_ari.csv", cfg), format = "file"),
  tar_target(sweep, run_sweep(proteins, meta, gene_sets, gmap, cfg)),
  tar_target(sweep_csv, write_out(sweep, "results", "sweep.csv", cfg), format = "file"),
  tar_target(fig_miss, save_fig(fig_missingness(qc_miss), "qc_missingness", cfg), format = "file"),
  tar_target(fig_cv, save_fig(fig_cv_performance(sweep), "cv_performance", cfg), format = "file"),
  tar_target(session_txt, write_session_info(cfg), format = "file")
)

if (quarto_ready) {
  targets <- c(targets, list(tar_quarto(report, path = here::here("t0", "report.qmd"))))
}

targets
