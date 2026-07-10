# HRvLR_F06_setup.R - Shared setup for Figure 6 (WGCNA module atlas).
# Unsupervised WGCNA on the FULL imputed proteome, then two module-level read-outs:
#   - baseline (T1) prediction of each training-adaptation trait (LOO cross-validated)
#   - responder (HR vs LR) signal per timepoint.
# Provides: wgcna_mem, wgcna_eig, pred, resp, pheno_tbl, pi_set, RPT_DIR, DAT_DIR.
# Plus all style.R / pathway_utils.R exports.

pacman::p_load(here, tidyverse, patchwork, grid)

source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))
source(here("04_Figures", "F06", "a_script", "functions", "clustering.R"))

RPT_DIR <- here("04_Figures", "F06", "b_reports")
DAT_DIR <- here("04_Figures", "F06", "c_data")
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

ci <- load_clustering_inputs()
pi_set <- ci$pi_set
pheno_tbl <- build_phenotype_table(here("00_input", "HRvLR_meta.csv"))

# Display labels for the phenotype traits (shared by panels C and S)
trait_lab <- c(
  comp_hypertrophy = "Comp. hypertrophy", d_mcsa = "d mCSA", d_fcsa_I = "d fCSA I",
  d_fcsa_II = "d fCSA II", d_1rm_legpress = "d 1RM leg", d_1rm_ext = "d 1RM ext"
)

wg <- run_wgcna(ci$full_abund, ci$meta)
wgcna_mem <- wg$membership
wgcna_eig <- wg$eigengene
# Network internals kept in memory for the WGCNA-QC supplements (soft-threshold
# sweep, gene dendrogram); not persisted, so c_data stays limited to the tables.
wgcna_sft <- wg$sft
wgcna_net <- wg$net
wgcna_power <- wg$chosen_power
wgcna_colors <- wg$module_colors

pred <- run_module_prediction(wgcna_eig, pheno_tbl)
resp <- run_module_responder(wgcna_eig)
pheno_lmm <- run_module_phenotype_lmm(wgcna_eig, pheno_tbl)

# WGCNA/AnnotationDbi mask tidyverse verbs; re-attach dplyr/tidyr last so the
# panels' unqualified verbs resolve correctly.
suppressPackageStartupMessages({
  for (p in c("tidyr", "dplyr")) {
    pk <- paste0("package:", p)
    if (pk %in% search()) suppressWarnings(detach(pk, character.only = TRUE))
    library(p, character.only = TRUE)
  }
})

write.csv(wgcna_mem, file.path(DAT_DIR, "wgcna_membership.csv"), row.names = FALSE)
write.csv(wgcna_eig, file.path(DAT_DIR, "wgcna_eigengene.csv"), row.names = FALSE)
write.csv(pheno_tbl, file.path(DAT_DIR, "phenotype.csv"), row.names = FALSE)
write.csv(pred, file.path(DAT_DIR, "module_prediction.csv"), row.names = FALSE)
write.csv(resp, file.path(DAT_DIR, "module_responder.csv"), row.names = FALSE)
write.csv(pheno_lmm, file.path(DAT_DIR, "module_phenotype_lmm.csv"), row.names = FALSE)

cat(sprintf(
  "Loaded: %d WGCNA modules (full proteome); prediction grid %d, responder grid %d\n",
  length(unique(wgcna_mem$group_id[wgcna_mem$group_id != "grey"])),
  nrow(pred), nrow(resp)
))
