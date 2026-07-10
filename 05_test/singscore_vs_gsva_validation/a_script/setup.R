# singscore_vs_gsva setup: score the same imputed matrix against the same
# Hallmark + GO Slim collection two ways. singscore ranks each sample on its own
# (single-sample, leakage-free); GSVA ranks each gene against the cohort
# (cohort-relative). Both matrices are trimmed to the pathways scored by both
# methods so the agreement analysis is like-for-like.

pacman::p_load(
  here, dplyr, readr, tibble, tidyr, ggplot2, limma, singscore, GSVA
)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))
set.seed(42)

rpt <- here("05_test", "singscore_vs_gsva_validation", "b_reports")
dat <- here("05_test", "singscore_vs_gsva_validation", "c_data")
dir.create(rpt, recursive = TRUE, showWarnings = FALSE)
dir.create(dat, recursive = TRUE, showWarnings = FALSE)

da <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
dep <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)

mat <- as.matrix(da$data)
sym <- dep$gene[match(rownames(mat), dep$uniprot_id)]
mat <- mat[!is.na(sym), ]
expr <- avereps(mat, ID = sym[!is.na(sym)])

pw_all <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
)
gene_sets <- pw_all[
  classify_database(names(pw_all)) %in% c("Hallmark", "GO Slim")
]

rank <- singscore::rankGenes(expr)
score_one <- function(genes) {
  present <- intersect(genes, rownames(expr))
  if (length(present) < 5) {
    return(rep(NA_real_, ncol(expr)))
  }
  singscore::simpleScore(rank, upSet = present)$TotalScore
}
sing <- t(vapply(gene_sets, score_one, numeric(ncol(expr))))
colnames(sing) <- colnames(expr)
sing <- sing[stats::complete.cases(sing), ]

gsva_par <- GSVA::gsvaParam(expr, gene_sets, minSize = 15, maxSize = 500)
gsva <- GSVA::gsva(gsva_par, verbose = FALSE)

common <- intersect(rownames(sing), rownames(gsva))
sing <- sing[common, colnames(expr)]
gsva <- gsva[common, colnames(expr)]

write_csv(
  as.data.frame(sing) |> rownames_to_column("pathway"),
  file.path(dat, "singscore_matrix.csv")
)
write_csv(
  as.data.frame(gsva) |> rownames_to_column("pathway"),
  file.path(dat, "gsva_recomputed.csv")
)

message(sprintf(
  "validation setup: %d pathways scored by both methods x %d samples",
  length(common), ncol(sing)
))
