# A feature bundle that has never been imputed.
#
# Every F06 cell reads the missForest matrix, which was imputed once across all
# 45 samples. A held-out subject's filled values were therefore informed by the
# subjects used to train the model that predicts them, which is information
# crossing the fold boundary before cross-validation starts. It is the channel
# that already sank the module -> d_mcsa result (04_Features/modules/loso_refit).
#
# The cheap test avoids imputation entirely rather than moving it in-fold: keep
# only proteins observed in every sample, and score those. 931 of 1900 qualify,
# so half the proteome survives the restriction. A lead that holds here was not
# built on filled values.
#
# Modules are absent by construction - an eigengene needs the WGCNA network, and
# that network is itself cohort-defined, which is the separate leak loso_refit
# measures.

pacman::p_load(here, dplyr, readr, limma)
source(here("functions", "shared_singscore.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "pred_features.R"))

complete_case_matrix <- function() {
  df <- read_csv(here("02_Normalization", "c_data", "normalized.csv"),
    show_col_types = FALSE
  )
  ann <- c("uniprot_id", "protein", "gene", "description")
  mat <- as.matrix(df[, setdiff(names(df), ann)])
  rownames(mat) <- df$uniprot_id
  keep <- rowSums(is.na(mat)) == 0
  list(mat = mat[keep, , drop = FALSE], gene = df$gene[keep])
}

pred_load_complete_case <- function() {
  cc <- complete_case_matrix()
  named <- !is.na(cc$gene) & cc$gene != ""
  expr <- avereps(cc$mat[named, , drop = FALSE], ID = cc$gene[named])

  pw <- build_pathway_collection(
    min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
  )
  gene_sets <- pw[classify_database(names(pw)) %in% c("Hallmark", "GO Slim")]

  da <- readRDS(pred_paths()$dalist)
  meta <- da$metadata |>
    transmute(
      sample = .data$Col_ID, subject = .data$Subject_ID,
      group = factor(.data$Group, levels = c("LR", "HR")),
      timepoint = factor(.data$Timepoint, levels = c("T1", "T2", "T3"))
    )

  list(
    meta = meta,
    feature_sets = list(
      proteins = expr,
      singscore = suppressWarnings(score_singscore(expr, gene_sets, min_size = 1L))
    ),
    pheno = read_csv(pred_paths()$pheno, show_col_types = FALSE)
  )
}
