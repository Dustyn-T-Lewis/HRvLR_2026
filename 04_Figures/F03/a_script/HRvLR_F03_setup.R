# HRvLR_F03_setup.R — Shared setup for Figure 3 (Clustering)
# Provides: wgcna_mem, wgcna_eig, mfuzz_mem, mfuzz_eig, spls_cv,
#           ora_results, concordance, pheno_models, perm_null, pheno_tbl,
#           gap_mat, pi_set, RPT_DIR, DAT_DIR
# Plus all style.R / pathway_utils.R exports (palettes, themes, helpers)

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(grid)
})

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

# Paths
WGCNA_MEM_FILE <- "05_Clustering/a_wgcna_paired/c_data/membership.csv"
WGCNA_EIG_FILE <- "05_Clustering/a_wgcna_paired/c_data/eigengene.csv"
MFUZZ_MEM_FILE <- "05_Clustering/b_mfuzz_gap/c_data/membership.csv"
MFUZZ_EIG_FILE <- "05_Clustering/b_mfuzz_gap/c_data/eigengene.csv"
SPLS_CV_FILE <- "05_Clustering/c_supervised/c_data/cv_error.csv"
ORA_FILE <- "05_Clustering/d_integration/c_data/ora_results.csv"
CONCORDANCE_FILE <- "05_Clustering/d_integration/c_data/concordance.csv"
PHENO_MOD_FILE <- "05_Clustering/d_integration/c_data/phenotype_models.csv"
PERM_NULL_FILE <- "05_Clustering/d_integration/c_data/perm_null.csv"
PHENO_FILE <- "05_Clustering/d_integration/c_data/phenotype.csv"

RPT_DIR <- "04_Figures/F03/b_reports"
DAT_DIR <- "04_Figures/F03/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load clustering outputs
wgcna_mem <- read.csv(WGCNA_MEM_FILE)
wgcna_eig <- read.csv(WGCNA_EIG_FILE)
mfuzz_mem <- read.csv(MFUZZ_MEM_FILE)
mfuzz_eig <- read.csv(MFUZZ_EIG_FILE)
spls_cv <- read.csv(SPLS_CV_FILE)
ora_results <- read.csv(ORA_FILE)
concordance <- read.csv(CONCORDANCE_FILE)
pheno_models <- read.csv(PHENO_MOD_FILE)
perm_null <- read.csv(PERM_NULL_FILE)
pheno_tbl <- read.csv(PHENO_FILE)

# HR-LR gap matrix + pi-gated set from shared clustering inputs
suppressMessages({
  source("05_Clustering/_shared_inputs.R")
  .ci <- load_clustering_inputs()
  gap_mat <- .ci$gap
  pi_set <- .ci$pi_set
})

# Timepoint order helper
tp_levels <- function() c("T1", "T2", "T3")

cat(sprintf(
  "Loaded: %d WGCNA modules, %d Mfuzz clusters, gap %dx%d, pi_set %d\n",
  length(unique(wgcna_mem$group_id[wgcna_mem$group_id != "grey"])),
  length(unique(mfuzz_mem$group_id)),
  nrow(gap_mat), ncol(gap_mat),
  length(pi_set)
))
