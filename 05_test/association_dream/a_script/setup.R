# association_dream setup: shared inputs for the dream mixed-model fits.
# Collapses the imputed protein matrix to gene symbols, scores each sample against
# a Hallmark + GO Slim gene-set collection with singscore, and assembles the design
# metadata (subject, group, timepoint), the subject-level cumulative volume load,
# and the continuous response phenotypes. Provides expr, scores, meta, gene_sets,
# cont_traits, and the report/data directories.

pacman::p_load(
  here, dplyr, readr, tibble, tidyr, ggplot2, limma, singscore, GSEABase
)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))
set.seed(42)

rpt <- here("05_test", "association_dream", "b_reports")
dat <- here("05_test", "association_dream", "c_data")
dir.create(rpt, recursive = TRUE, showWarnings = FALSE)
dir.create(dat, recursive = TRUE, showWarnings = FALSE)

da <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
dep <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)
pheno <- read_csv(
  here("04_Figures", "F06", "c_data", "phenotype.csv"),
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
scores <- t(vapply(gene_sets, score_one, numeric(ncol(expr))))
colnames(scores) <- colnames(expr)
scores <- scores[stats::complete.cases(scores), ]

write_csv(
  as.data.frame(scores) |> rownames_to_column("pathway"),
  file.path(dat, "singscore_matrix.csv")
)

accum_vl <- da$metadata |>
  filter(!is.na(ACCUM_VL)) |>
  distinct(Subject_ID, ACCUM_VL)

cont_traits <- c(
  comp_hypertrophy = "Composite hypertrophy",
  d_fcsa_I = "delta fCSA type I",
  d_fcsa_II = "delta fCSA type II",
  d_mcsa = "delta mCSA",
  d_1rm_legpress = "delta 1RM leg press",
  d_1rm_ext = "delta 1RM ext"
)

meta <- da$metadata |>
  transmute(
    sample = Col_ID, subject = Subject_ID,
    group = factor(Group, levels = c("LR", "HR")),
    timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  ) |>
  left_join(accum_vl, by = c("subject" = "Subject_ID")) |>
  left_join(
    pheno |> dplyr::select(subject, all_of(names(cont_traits))),
    by = "subject"
  ) |>
  mutate(accum_vl_z = as.numeric(scale(ACCUM_VL)))
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample
meta <- meta[colnames(scores), ]

message(sprintf(
  "association_dream setup: %d genes, %d pathways, %d samples, %d VL missing",
  nrow(expr), nrow(scores), ncol(scores), sum(is.na(meta$ACCUM_VL))
))
