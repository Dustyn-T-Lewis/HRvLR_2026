# 05_test setup - GSVA pathway scores + mixed models (experimental).
# Turns the imputed protein matrix into per-sample gene-symbol expression, loads
# the Hallmark + GO Slim gene sets, the design metadata, and the response
# phenotypes. Provides expr, GENE_SETS, meta, pheno, dirs.

pacman::p_load(here, dplyr, readr, tibble, tidyr, ggplot2, limma, GSVA, lme4, lmerTest)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))
set.seed(42)

RPT <- here("05_test", "b_reports")
DAT <- here("05_test", "c_data")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

da <- readRDS(here("02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"))
dep <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)
pheno <- read_csv(here("04_Figures", "F06", "c_data", "phenotype.csv"), show_col_types = FALSE)

meta <- da$metadata |>
  transmute(
    sample = Col_ID, subject = Subject_ID,
    group = factor(Group, levels = c("LR", "HR")),
    timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  )

mat <- as.matrix(da$data)
sym <- dep$gene[match(rownames(mat), dep$uniprot_id)]
mat <- mat[!is.na(sym), ]
expr <- avereps(mat, ID = sym[!is.na(sym)])

pw_all <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
)
GENE_SETS <- pw_all[classify_database(names(pw_all)) %in% c("Hallmark", "GO Slim")]

cat(sprintf(
  "05_test setup: %d genes x %d samples, %d gene sets\n",
  nrow(expr), ncol(expr), length(GENE_SETS)
))
