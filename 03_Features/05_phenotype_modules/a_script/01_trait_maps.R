#!/usr/bin/env Rscript
# Label-free association maps: do the WGCNA eigengenes or the singscore
# pathway scores track the continuous phenotypes, tested independently at
# each timepoint?
#
# The HR/LR label appears nowhere. Each feature (12 module eigengenes,
# 57 Hallmark singscores) is correlated with each of the six phenotype
# outcomes at T1, T2 and T3 separately, Spearman against a B = 1000
# subject-permutation null with one shared permutation index, so
# replicate b is the same shuffled cohort everywhere. A feature counts
# consistent for a phenotype when the sign agrees at all three
# timepoints and empirical p < 0.05 at two or more; the same call runs
# inside every permuted cohort, giving the consistency count its own
# null. Timepoints share subjects, so consistency is stability, not
# independent replication — and the phenotypes are training changes, so
# T1 associations are baseline forecasts, which F06 already showed do
# not exist at the protein level.
#
# Read 01_{module,pathway}_trait.csv for the full maps and
# 01_consistency.csv for observed-vs-null consistency counts per
# phenotype. A null looks like observed counts inside their permuted
# distribution, which is what n = 16 should produce.

pacman::p_load(here, dplyr, tidyr, readr, tibble, purrr, cluster)
source(here(
  "03_Features", "05_phenotype_modules", "a_script", "trait_helpers.R"
))

set.seed(42)
OUT <- here("03_Features", "05_phenotype_modules", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TIMEPOINTS <- c("T1", "T2", "T3")
PHENOS <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

pheno <- read_csv(
  here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
) |>
  mutate(subject = sub("^(HR|LR)_", "", .data$subject))

split_id <- function(ids) {
  tibble(
    subject = sub("_T\\d$", "", ids),
    timepoint = sub("^.*_", "", ids)
  )
}

eig <- read_csv(
  here("03_Features", "03_WGCNA", "c_data", "wgcna_eigengene.csv"),
  show_col_types = FALSE
)
eig <- bind_cols(eig, rename(split_id(eig$sample_id), subj = "subject"))

module_mats <- map(TIMEPOINTS, function(tp) {
  eig |>
    filter(.data$timepoint == tp) |>
    select("group_id", "subj", "ME") |>
    pivot_wider(names_from = "subj", values_from = "ME") |>
    column_to_rownames("group_id") |>
    as.matrix()
})

sing <- readRDS(
  here("03_Features", "02_Pathways", "c_data", "singscore_scores.rds")
)
sing_meta <- split_id(colnames(sing))
pathway_mats <- map(TIMEPOINTS, function(tp) {
  m <- sing[, sing_meta$timepoint == tp, drop = FALSE]
  colnames(m) <- sing_meta$subject[sing_meta$timepoint == tp]
  m
})

complete_subj <- reduce(map(module_mats, colnames), intersect)
message(sprintf(
  "%d of %d subjects have all three timepoints", length(complete_subj),
  nrow(pheno)
))
module_mats <- map(module_mats, \(m) m[, complete_subj, drop = FALSE])
pathway_mats <- map(pathway_mats, \(m) m[, complete_subj, drop = FALSE])
pheno <- pheno[match(complete_subj, pheno$subject), ]

perm_idx <- perm_index(length(complete_subj))

scan_level <- function(mats, level) {
  per_pheno <- map(PHENOS, function(ph) {
    scans <- map(mats, cor_scan, pheno = pheno[[ph]], perm_idx = perm_idx)
    cons <- consistency_scan(scans)
    list(
      map = imap_dfr(scans, \(s, i) mutate(
        s$obs,
        timepoint = TIMEPOINTS[i], phenotype = ph
      )),
      consistency = tibble(
        level = level, phenotype = ph,
        n_observed = cons$n_observed,
        null_median = cons$null_median, emp_p = cons$emp_p,
        consistent = paste(cons$consistent, collapse = ";")
      )
    )
  })
  list(
    map = bind_rows(map(per_pheno, "map")) |> mutate(level = level),
    consistency = bind_rows(map(per_pheno, "consistency"))
  )
}

modules <- scan_level(module_mats, "module")
pathways <- scan_level(pathway_mats, "pathway")

write_csv(modules$map, file.path(OUT, "01_module_trait.csv"))
write_csv(pathways$map, file.path(OUT, "01_pathway_trait.csv"))
consistency <- bind_rows(modules$consistency, pathways$consistency)
write_csv(consistency, file.path(OUT, "01_consistency.csv"))

message(paste(
  capture.output(print(as.data.frame(
    select(consistency, -"consistent")
  ))),
  collapse = "\n"
))
