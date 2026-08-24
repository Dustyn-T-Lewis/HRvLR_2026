# Per-sample pathway scores, computed once and read by every downstream stage.
#
# This script exists because V1 did not have one. `singscore_scores.rds` was
# tracked there with no code that produced it: three files read it and none
# wrote it, the same orphan V1's own 00_input/a_script/01_build_phenotype.R
# header was added to fix for phenotype.csv. Scores that no script can rebuild
# cannot be checked, so V2 rebuilds them.
#
# singscore ranks within each sample, so a score depends only on that sample's
# own gene ranks and carries no information from the rest of the cohort. That is
# what makes it safe to compute here, upstream of any label, and reuse
# everywhere without recomputing per contrast.

pacman::p_load(here, dplyr)

source(here("functions", "shared_singscore.R"))
source(here("functions", "shared_pathway_utils.R"))

OUT <- here("02_Normalization", "c_data", "singscore_scores.rds")

imputed <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
expr <- as.matrix(imputed$data)
rownames(expr) <- imputed$annotation$gene
expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", , drop = FALSE]

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |>
  (\(d) split(d$gene_symbol, d$gs_name))() |>
  lapply(unique)
goslim <- build_goslim_gene_sets(min_size = 10, max_size = 500)

scores <- score_singscore(c(hallmark, goslim),
  expr = expr,
  min_size = SET_FLOOR
)

saveRDS(scores, OUT)

message(sprintf(
  "%d pathway sets x %d samples (%d Hallmark, %d GO Slim after the %d-member floor)",
  nrow(scores), ncol(scores),
  sum(grepl("^HALLMARK", rownames(scores))),
  sum(grepl("^GOSLIM", rownames(scores))), SET_FLOOR
))
