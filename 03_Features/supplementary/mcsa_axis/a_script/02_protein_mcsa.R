#!/usr/bin/env Rscript
# Does whole-muscle CSA change have a protein-level correlate of its own?
#
# The 931 complete-case proteins against d_mcsa, blood index partialled out,
# fitted separately at each timepoint and on the training delta. limma rather
# than a mixed model because a config holds one sample per subject, so there
# are no repeated measures left to block on; the moderated variance is what
# limma is here for at n = 14-15. BH within a config, never across the four,
# since the configs share subjects and are not independent tests.
#
# Only T1 precedes the outcome d_mcsa measures, so only T1 can be a forecast.
# Every config runs twice, with LR_S14 and without: that subject posts the
# study's largest whole-muscle gain alongside near-worst fibre change, and the
# two runs side by side are the result rather than one being the answer.
#
# Read 02_protein_mcsa.csv for the per-protein slopes and BH q by config and
# subject set, 02_survivors.csv for anything clearing BH, and 02_permutation.csv
# if the permutation armed. A null looks like no BH survivor in any config,
# which at n = 15 it should.

pacman::p_load(here, dplyr, purrr, readr, tibble)
source(here("functions", "shared_hlm.R"))
source(here("functions", "blood_index_model.R"))
source(here(
  "03_Features", "supplementary", "galamm_pilot", "a_script",
  "pilot_helpers.R"
))
source(here(
  "03_Features", "supplementary", "mcsa_axis", "a_script",
  "mcsa_helpers.R"
))

OUT <- here("03_Features", "supplementary", "mcsa_axis", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_PERM <- 200L

inp <- pilot_data()
pheno <- phenotype_table()

subject_sets <- list(
  `all subjects` = NULL, `LR_S14 dropped` = DISCORDANT_SUBJECT
)

scans <- expand.grid(
  config = MCSA_CONFIGS, subject_set = names(subject_sets),
  stringsAsFactors = FALSE
) |>
  pmap(function(config, subject_set) {
    design <- config_design(
      inp$mat, inp$meta, pheno, config,
      drop = subject_sets[[subject_set]]
    )
    mcsa_scan(design) |>
      mutate(
        config = config, subject_set = subject_set,
        n = length(design$subject)
      )
  }) |>
  list_rbind()

write_csv(scans, file.path(OUT, "02_protein_mcsa.csv"))

survivors <- filter(scans, .data$bh < 0.05)
write_csv(survivors, file.path(OUT, "02_survivors.csv"))

if (nrow(survivors)) {
  perm <- survivors |>
    group_by(.data$config, .data$subject_set) |>
    group_split() |>
    map(function(hits) {
      design <- config_design(
        inp$mat, inp$meta, pheno, hits$config[1],
        drop = subject_sets[[hits$subject_set[1]]]
      )
      mcsa_permutation(design, hits, n_perm = N_PERM) |>
        mutate(config = hits$config[1], subject_set = hits$subject_set[1])
    }) |>
    list_rbind()
  write_csv(perm, file.path(OUT, "02_permutation.csv"))

  checks <- survivors |>
    pmap(function(feature, config, subject_set, ...) {
      design <- config_design(
        inp$mat, inp$meta, pheno, config,
        drop = subject_sets[[subject_set]]
      )
      survivor_checks(design, feature) |>
        mutate(config = config, subject_set = subject_set)
    }) |>
    list_rbind()
  write_csv(checks, file.path(OUT, "02_survivor_checks.csv"))
  print(as.data.frame(checks), digits = 3, row.names = FALSE)
}

summary_by_cell <- scans |>
  summarise(
    n = first(.data$n), min_bh = min(.data$bh), nominal = sum(.data$p < 0.05),
    expected = 0.05 * n(), survivors = sum(.data$bh < 0.05),
    .by = c("config", "subject_set")
  )
write_csv(summary_by_cell, file.path(OUT, "02_scan_summary.csv"))

print(as.data.frame(summary_by_cell), digits = 3, row.names = FALSE)
cat(sprintf(
  "\n%d BH survivors across %d cells; permutation %s\n",
  nrow(survivors), nrow(summary_by_cell),
  if (nrow(survivors)) "armed" else "never armed"
))
