# F02 setup: the normalized, imputed, and DEP artifacts from upstream c_data, plus
# the metadata and the contrast names. Panels and run.R source this first. Writes nothing.
# Provides: norm_df, imp_df, dep_df, meta, samp_names (norm), imp_samps (imp),
#           ann_cols, imp_mat, MAIN_CONTRASTS, RPT_DIR, DAT_DIR
# Plus all style.R exports (palettes, themes, helpers)

pacman::p_load(here, tidyverse, patchwork, grid, vegan)

source(here("04_Figures", "functions", "style.R"))
source(here("03_DEP", "contrasts.R"))

# Paths
NORM_FILE <- here("02_Normalization", "c_data", "normalized.csv")
IMP_FILE <- here("02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds")
DEP_FILE <- here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv")
META_FILE <- here("00_input", "HRvLR_meta.csv")

if (!exists("F02_AUDIT")) F02_AUDIT <- list()

RPT_DIR <- here("04_Figures", "F02_proteome", "b_reports")
DAT_DIR <- here("04_Figures", "F02_proteome", "c_data")

# The between-responder-at-timepoint contrasts (Trained_HRvLR, Acute_HRvLR) are
# intentionally dropped; F02 reads responder divergence from the interaction terms.
# Names validated against the single source so the subset can't drift.
MAIN_CONTRASTS <- c(
  "Baseline_HRvLR",
  "Training_HR", "Training_LR", "Training_Interaction",
  "Acute_HR", "Acute_LR", "Acute_Interaction"
)
stopifnot(all(MAIN_CONTRASTS %in% sub(" =.*$", "", HRVLR_CONTRASTS)))

# Load data
norm_df <- read_csv(NORM_FILE, show_col_types = FALSE)
imp_dal <- readRDS(IMP_FILE)
imp_df <- bind_cols(
  as_tibble(imp_dal$annotation[, c("uniprot_id", "protein", "gene", "description")]),
  as_tibble(imp_dal$data)
)
dep_df <- read_csv(DEP_FILE, show_col_types = FALSE)

ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

# Imputed data may have fewer samples (QC-removed during imputation)
imp_ann <- intersect(names(imp_df), ann_cols)
imp_samps <- setdiff(names(imp_df), imp_ann)

# Metadata - from CSV, joined to imp_samps (the analysis-ready sample set)
meta_raw <- read_csv(META_FILE, show_col_types = FALSE)
meta <- tibble(sample_id = imp_samps) |>
  left_join(meta_raw |> select(Col_ID, Subject_ID, Group, Timepoint, Group_Time),
    by = c("sample_id" = "Col_ID")
  ) |>
  mutate(
    Col_ID = sample_id,
    subject = Subject_ID,
    group = factor(Group_Time,
      levels = c(
        "HR_T1", "HR_T2", "HR_T3",
        "LR_T1", "LR_T2", "LR_T3"
      )
    ),
    Group = factor(Group, levels = c("HR", "LR")),
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  )

cat(sprintf(
  "Loaded: %d norm proteins (%d samples), %d imp proteins (%d samples), %d DEP rows\n",
  nrow(norm_df), length(samp_names), nrow(imp_df), length(imp_samps), nrow(dep_df)
))

# Imputed matrix (proteins x samples)
imp_mat <- as.matrix(imp_df[, imp_samps])
rownames(imp_mat) <- imp_df$gene
