#!/usr/bin/env Rscript
# Do proteins that move together over training form clusters that track
# the phenotype outcomes?
#
# WGCNA clusters abundance co-expression; this clusters the training
# response itself. For every subject with both biopsies the T2 - T1
# delta is computed on the 931 complete-case proteins (nothing imputed),
# z-scored per protein, and clustered by hclust ward.D2 on correlation
# distance with k chosen by mean silhouette over 2:10 — all before any
# phenotype is read. Each cluster's per-subject mean delta is then
# correlated with the six phenotype outcomes, Spearman against the same
# B = 1000 subject-permutation null as 01_trait_maps.R.
#
# Read 02_clusters.csv for the protein memberships and
# 02_cluster_trait.csv for the association table. A null looks like no
# cluster-phenotype pair beating its permutation p, and with k chosen
# unsupervised the cluster count is a description, not a result.

pacman::p_load(here, dplyr, tidyr, readr, tibble, purrr, cluster)
source(here("functions", "shared_hlm.R"))
source(here("functions", "blood_index_model.R"))
source(here(
  "03_Features", "supplementary", "galamm_pilot", "a_script", "pilot_helpers.R"
))
source(here(
  "03_Features", "supplementary", "phenotype_modules", "a_script",
  "trait_helpers.R"
))

set.seed(42)
OUT <- here("03_Features", "supplementary", "phenotype_modules", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PHENOS <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

inp <- pilot_data()
meta <- inp$meta
meta$subj <- sub("^(HR|LR)_", "", as.character(meta$subject))

paired <- intersect(
  meta$subj[meta$timepoint == "T1"], meta$subj[meta$timepoint == "T2"]
)
delta <- vapply(paired, function(s) {
  inp$mat[, meta$sample[meta$subj == s & meta$timepoint == "T2"]] -
    inp$mat[, meta$sample[meta$subj == s & meta$timepoint == "T1"]]
}, numeric(nrow(inp$mat)))
message(sprintf(
  "%d subjects with paired T1/T2 biopsies", length(paired)
))

delta_z <- t(scale(t(delta)))
d <- as.dist(1 - cor(t(delta_z)))
hc <- hclust(d, method = "ward.D2")
k <- choose_k(d, hc)
membership <- cutree(hc, k)
message(sprintf("silhouette selects k = %d", k))

cluster_scores <- vapply(seq_len(k), function(cl) {
  colMeans(delta_z[membership == cl, , drop = FALSE])
}, numeric(ncol(delta_z)))
colnames(cluster_scores) <- paste0("C", seq_len(k))

pheno <- read_csv(
  here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
) |>
  mutate(subject = sub("^(HR|LR)_", "", .data$subject))
pheno <- pheno[match(paired, pheno$subject), ]

perm_idx <- perm_index(length(paired))
cluster_trait <- map_dfr(PHENOS, function(ph) {
  cor_scan(t(cluster_scores), pheno[[ph]], perm_idx)$obs |>
    mutate(phenotype = ph)
})

anno <- readRDS(
  here("02_Normalization", "c_data", "DAList_normalized.rds")
)$annotation
write_csv(
  tibble(
    feature = rownames(delta_z),
    gene = anno$gene[match(rownames(delta_z), anno$uniprot_id)],
    cluster = paste0("C", membership)
  ),
  file.path(OUT, "02_clusters.csv")
)
write_csv(cluster_trait, file.path(OUT, "02_cluster_trait.csv"))

message(paste(
  capture.output(print(as.data.frame(
    filter(cluster_trait, .data$emp_p < 0.05)
  ))),
  collapse = "\n"
))
