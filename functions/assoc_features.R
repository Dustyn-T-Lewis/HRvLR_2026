# Shared inputs and feature builders for the F04 association suite. Association
# is in-sample and descriptive, so it reads the normalized (not imputed)
# proteome and keeps NA; limma and the mixed model handle missing values per
# feature. Three feature spaces: proteins (uniprot-level, gene-labelled to match
# 03_DEP), pathway scores (singscore), and WGCNA module eigengenes.

pacman::p_load(here, dplyr, tibble, readr, limma)
source(here("functions", "shared_hlm.R"))
source(here("functions", "pred_features.R"))
source(here("04_Figures", "functions", "f04_association.R"))

assoc_paths <- function() {
  list(
    norm = here("02_Normalization", "c_data", "DAList_normalized.rds"),
    imputed = here(
      "02_Normalization", "imputation", "c_data",
      "DAList_imputed_missforest.rds"
    ),
    pheno = here("00_input", "c_data", "phenotype.csv"),
    eigen = here(
      "04_Figures", "F04_association", "WGCNA", "c_data", "wgcna_eigengene.csv"
    ),
    dep = here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv")
  )
}

# Load once: sample metadata, the three association feature matrices (features x
# sample), the gene label map, the phenotype table, and an imputed matrix
# for the variance partition (which needs an NA-free input).
assoc_load <- function() {
  paths <- assoc_paths()
  da <- readRDS(paths$norm)
  meta <- hlm_meta(da)

  proteins <- as.matrix(da$data)
  gene_map <- setNames(da$annotation$gene, da$annotation$uniprot_id)
  imp <- readRDS(paths$imputed)

  # Pathways reuse the exact leakage-free singscore matrix the F05 suite scores
  # (Hallmark + GO Slim on the completed proteome), so the two figures share one
  # pathway space. Proteins stay on the measured normalized values for the DE.
  pathways <- pred_singscore(
    pred_gene_expression(imp), pred_paths()$cache
  )[, colnames(proteins)]

  feature_sets <- list(
    proteins = proteins,
    pathways = pathways,
    modules = pred_eigengene_matrix(paths$eigen, colnames(proteins))
  )

  list(
    meta = meta, meta_da = da$metadata, feature_sets = feature_sets,
    gene_map = gene_map, pheno = read_csv(paths$pheno, show_col_types = FALSE),
    imputed_proteins = as.matrix(imp$data)
  )
}

# Leaner per-contrast group difference: moderated limma `~ group` on the
# contrast's per-subject matrix (levels for T1/T2/T3, changes for training and
# acute). The snapshot view that branch 02 shows beside the mixed-model slice.
group_limma_contrast <- function(feature_mat, meta, contrast, group_of) {
  x <- t(pred_contrast_matrix(feature_mat, meta, contrast))
  g <- factor(group_of[colnames(x)], levels = c("LR", "HR"))
  fit <- eBayes(lmFit(x, model.matrix(~g)))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tibble(
    feature = rownames(tt), effect = tt$logFC, t = tt$t,
    p = tt$P.Value, bh = tt$adj.P.Val, n = ncol(x)
  )
}

# HR/LR label per subject, from the sample metadata.
assoc_group_of <- function(meta) {
  g <- meta |> dplyr::distinct(subject, group)
  setNames(as.character(g$group), g$subject)
}

# 03_DEP contrast name for each branch-02 contrast, so the reconciliation reads
# the one source of truth rather than recomputing it.
DEP_CONTRAST <- c(
  T1 = "Baseline_HRvLR", T2 = "Trained_HRvLR", T3 = "Acute_HRvLR",
  training = "Training_Interaction", acute = "Acute_Interaction"
)

dep_reference <- function(contrast) {
  paths <- assoc_paths()
  dep <- readr::read_csv(paths$dep, show_col_types = FALSE)
  nm <- DEP_CONTRAST[[contrast]]
  tibble(
    feature = dep$uniprot_id,
    dep_logfc = dep[[paste0("logFC_", nm)]],
    dep_p = dep[[paste0("P.Value_", nm)]],
    dep_bh = dep[[paste0("adj.P.Val_", nm)]]
  )
}
