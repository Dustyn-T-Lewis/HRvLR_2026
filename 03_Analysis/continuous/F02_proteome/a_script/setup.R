# F02 setup, continuous tree: the normalized, imputed, and DEP artifacts,
# plus the phenotype table (for the colour gradient PCA uses in place of a
# group). Mirrors categorical/F02_proteome/a_script/setup.R with the group
# columns dropped and the two pooled contrasts in place of the nine.
# Panels and 01_run_proteome.R source this first. Writes nothing.

pacman::p_load(here, tidyverse, patchwork, grid, vegan)

source(here("functions", "shared_style.R"))
source(here("03_Features", "contrasts.R"))

NORM_FILE <- here("02_Normalization", "c_data", "normalized.csv")
IMP_FILE <- here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
)
DEP_FILE <- here(
  "03_Analysis", "continuous", "01_Proteins", "c_data",
  "02_combined_results.csv"
)
META_FILE <- here("00_input", "HRvLR_meta.csv")
PHENO_FILE <- here("00_input", "c_data", "phenotype.csv")

if (!exists("F02_AUDIT")) F02_AUDIT <- list()

RPT_DIR <- here("03_Analysis", "continuous", "F02_proteome", "b_reports")
DAT_DIR <- here("03_Analysis", "continuous", "F02_proteome", "c_data")

MAIN_CONTRASTS <- c("Training", "Acute")
stopifnot(all(MAIN_CONTRASTS %in% trimws(sub("=.*$", "", POOLED_CONTRASTS))))

norm_df <- read_csv(NORM_FILE, show_col_types = FALSE)
imp_dal <- readRDS(IMP_FILE)
imp_df <- bind_cols(
  as_tibble(
    imp_dal$annotation[, c("uniprot_id", "protein", "gene", "description")]
  ),
  as_tibble(imp_dal$data)
)
dep_df <- read_csv(DEP_FILE, show_col_types = FALSE)

ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

imp_ann <- intersect(names(imp_df), ann_cols)
imp_samps <- setdiff(names(imp_df), imp_ann)

meta_raw <- read_csv(META_FILE, show_col_types = FALSE)
pheno <- read_csv(PHENO_FILE, show_col_types = FALSE)

meta <- tibble(sample_id = imp_samps) |>
  left_join(
    meta_raw |> select(Col_ID, Subject_ID, Timepoint),
    by = c("sample_id" = "Col_ID")
  ) |>
  mutate(
    Col_ID = sample_id, subject = Subject_ID,
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  ) |>
  left_join(
    pheno |> select(subject, comp_hypertrophy),
    by = "subject"
  )

cat(sprintf(
  paste0(
    "Loaded: %d norm proteins (%d samples), %d imp proteins (%d samples), ",
    "%d DEP rows\n"
  ),
  nrow(norm_df), length(samp_names), nrow(imp_df), length(imp_samps),
  nrow(dep_df)
))

imp_mat <- as.matrix(imp_df[, imp_samps])
rownames(imp_mat) <- imp_df$gene
