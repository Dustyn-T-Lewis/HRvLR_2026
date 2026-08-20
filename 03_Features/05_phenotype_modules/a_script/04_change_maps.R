#!/usr/bin/env Rscript
# The change configurations the per-timepoint maps skipped: does the
# training response (T2 - T1) or the acute response (T3 - T2) of a
# module eigengene or pathway score track the phenotype outcomes?
#
# Same machinery as 01_trait_maps.R — Spearman against B = 1000
# subject-permutation nulls — applied to per-subject deltas instead of
# static timepoints. No design contrast is needed for this: the deltas
# are computed directly per subject, and only subjects with both
# biopsies of a config enter it. The two configs answer different
# physiology, so there is no cross-config consistency call; each map
# stands alone. The acute protein-delta clustering mirrors
# 02_delta_clusters.R for the T3 - T2 response.
#
# Read 04_change_trait.csv for the module and pathway maps and
# 04_acute_cluster{s,_trait}.csv for the acute clusters. A null looks
# like no feature beating its permutation p, which the per-timepoint
# maps already made the expected outcome.

pacman::p_load(here, dplyr, tidyr, readr, tibble, purrr, cluster)
source(here("functions", "shared_hlm.R"))
source(here("functions", "blood_index_model.R"))
source(here(
  "03_Features", "04_galamm_pilot", "a_script", "pilot_helpers.R"
))
source(here(
  "03_Features", "05_phenotype_modules", "a_script", "trait_helpers.R"
))

set.seed(42)
OUT <- here("03_Features", "05_phenotype_modules", "c_data")

PHENOS <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)
CONFIGS <- list(training = c("T1", "T2"), acute = c("T2", "T3"))

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

# features x samples matrix -> features x subjects delta matrix for one
# config, keeping subjects with both biopsies.
delta_mat <- function(mat, ids, config) {
  from <- config[1]
  to <- config[2]
  have <- split(ids$subject, ids$timepoint)
  subjects <- intersect(have[[from]], have[[to]])
  col_of <- function(tp) {
    match(paste(subjects, tp, sep = "_"), paste(
      ids$subject, ids$timepoint,
      sep = "_"
    ))
  }
  out <- mat[, col_of(to), drop = FALSE] - mat[, col_of(from), drop = FALSE]
  colnames(out) <- subjects
  out
}

scan_deltas <- function(mat, ids, level) {
  imap_dfr(CONFIGS, function(config, config_name) {
    dm <- delta_mat(mat, ids, config)
    ph <- pheno[match(colnames(dm), pheno$subject), ]
    perm_idx <- perm_index(ncol(dm))
    map_dfr(PHENOS, function(p) {
      cor_scan(dm, ph[[p]], perm_idx)$obs |>
        mutate(phenotype = p, config = config_name)
    })
  }) |>
    mutate(level = level)
}

eig <- read_csv(
  here(
    "03_Analysis", "categorical", "03_WGCNA", "c_data", "wgcna_eigengene.csv"
  ),
  show_col_types = FALSE
) |>
  pivot_wider(names_from = "sample_id", values_from = "ME")
eig_mat <- as.matrix(column_to_rownames(eig, "group_id"))

sing <- readRDS(here("02_Normalization", "c_data", "singscore_scores.rds"))

change_trait <- bind_rows(
  scan_deltas(eig_mat, split_id(colnames(eig_mat)), "module"),
  scan_deltas(sing, split_id(colnames(sing)), "pathway")
)
write_csv(change_trait, file.path(OUT, "04_change_trait.csv"))

inp <- pilot_data()
prot_ids <- tibble(
  subject = sub("^(HR|LR)_", "", as.character(inp$meta$subject)),
  timepoint = as.character(inp$meta$timepoint)
)
acute_delta <- delta_mat(inp$mat, prot_ids, CONFIGS$acute)
acute_z <- t(scale(t(acute_delta)))
d <- as.dist(1 - cor(t(acute_z)))
hc <- hclust(d, method = "ward.D2")
k <- choose_k(d, hc)
membership <- cutree(hc, k)
message(sprintf(
  "acute deltas: %d subjects, silhouette selects k = %d",
  ncol(acute_delta), k
))

cluster_scores <- vapply(seq_len(k), function(cl) {
  colMeans(acute_z[membership == cl, , drop = FALSE])
}, numeric(ncol(acute_z)))
colnames(cluster_scores) <- paste0("A", seq_len(k))

ph <- pheno[match(colnames(acute_delta), pheno$subject), ]
perm_idx <- perm_index(ncol(acute_delta))
acute_trait <- map_dfr(PHENOS, function(p) {
  cor_scan(t(cluster_scores), ph[[p]], perm_idx)$obs |>
    mutate(phenotype = p)
})

anno <- readRDS(
  here("02_Normalization", "c_data", "DAList_normalized.rds")
)$annotation
write_csv(
  tibble(
    feature = rownames(acute_z),
    gene = anno$gene[match(rownames(acute_z), anno$uniprot_id)],
    cluster = paste0("A", membership)
  ),
  file.path(OUT, "04_acute_clusters.csv")
)
write_csv(acute_trait, file.path(OUT, "04_acute_cluster_trait.csv"))

hits <- change_trait |>
  filter(.data$emp_p < 0.05) |>
  count(.data$level, .data$config, .data$phenotype)
message(paste(
  capture.output(print(as.data.frame(hits))),
  collapse = "\n"
))
message(paste(
  capture.output(print(as.data.frame(
    filter(acute_trait, .data$emp_p < 0.05)
  ))),
  collapse = "\n"
))
