# Shared inputs and feature/outcome builders for the F05 prediction suite.
# Loads the imputed proteome, exposes three feature spaces (gene-collapsed
# proteins, singscore pathway scores, WGCNA module eigengenes), and builds the
# per-contrast subject-level feature matrices plus the phenotype outcomes.
#
# singscore is single-sample and rank-based: each sample's score depends only on
# its own within-sample gene ranks, so scoring the full matrix once is
# leakage-free. Proteins are per-sample values, also leakage-free. Module
# eigengenes are cohort-relative (fit across all samples), so they are flagged
# on the figure as an optimistic upper bound rather than refit per fold. Every
# fold-specific transform (z-scoring, tuning) still happens train-only in the
# harness.

pacman::p_load(here, dplyr, readr, limma, singscore)
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_singscore.R"))

pred_paths <- function() {
  list(
    dalist = here(
      "02_Normalization", "imputation", "c_data",
      "DAList_imputed_missforest.rds"
    ),
    pheno = here("00_input", "c_data", "phenotype.csv"),
    eigen = here(
      "04_Figures", "F04_association", "WGCNA", "c_data", "wgcna_eigengene.csv"
    ),
    cache = here(
      "04_Figures", "F05_prediction", "shared", "c_data",
      "singscore_scores.rds"
    )
  )
}

pred_gene_expression <- function(da) {
  mat <- as.matrix(da$data)
  sym <- da$annotation$gene[match(rownames(mat), da$annotation$uniprot_id)]
  keep <- !is.na(sym) & sym != ""
  avereps(mat[keep, , drop = FALSE], ID = sym[keep])
}

pred_singscore <- function(expr, cache_path) {
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  pw <- build_pathway_collection(
    min_size = 15, max_size = 500,
    include_goslim = TRUE, exclude_variants = TRUE
  )
  gene_sets <- pw[classify_database(names(pw)) %in% c("Hallmark", "GO Slim")]
  scores <- suppressWarnings(score_singscore(expr, gene_sets, min_size = 1L))
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(scores, cache_path)
  scores
}

pred_eigengene_matrix <- function(eigen_path, sample_order) {
  long <- read_csv(eigen_path, show_col_types = FALSE)
  mods <- sort(unique(long$group_id))
  wide <- matrix(NA_real_,
    nrow = length(mods), ncol = length(sample_order),
    dimnames = list(paste0("ME_", mods), sample_order)
  )
  wide[cbind(paste0("ME_", long$group_id), long$sample_id)] <- long$ME
  wide
}

# Load everything once and return the shared bundle: sample metadata, the three
# feature matrices (features x sample), and the phenotype table.
pred_load <- function() {
  paths <- pred_paths()
  da <- readRDS(paths$dalist)

  meta <- da$metadata |>
    transmute(
      sample = Col_ID, subject = Subject_ID,
      group = factor(Group, levels = c("LR", "HR")),
      timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
    )

  expr <- pred_gene_expression(da)
  feature_sets <- list(
    proteins   = expr,
    singscore  = pred_singscore(expr, paths$cache),
    eigengenes = pred_eigengene_matrix(paths$eigen, colnames(da$data))
  )

  list(
    meta = meta, feature_sets = feature_sets,
    pheno = read_csv(paths$pheno, show_col_types = FALSE)
  )
}

pred_sample_at <- function(meta, subj, tp) {
  hit <- meta$sample[meta$subject == subj & meta$timepoint == tp]
  if (length(hit) == 1) hit else NA_character_
}

# Subject-level feature matrix (subjects x features) for one contrast. T1/T2/T3
# and baseline (= T1) are the level at that timepoint; training is the T2 - T1
# change, acute the T3 - T2 change. Subjects missing a required sample drop out.
pred_contrast_matrix <- function(feature_mat, meta, contrast) {
  spec <- switch(contrast,
    baseline = list(tps = "T1", diff = FALSE),
    T1 = list(tps = "T1", diff = FALSE),
    T2 = list(tps = "T2", diff = FALSE),
    T3 = list(tps = "T3", diff = FALSE),
    training = list(tps = c("T1", "T2"), diff = TRUE),
    acute = list(tps = c("T2", "T3"), diff = TRUE),
    stop("unknown contrast: ", contrast)
  )
  subjects <- levels(factor(meta$subject))
  rows <- lapply(subjects, function(s) {
    ids <- vapply(
      spec$tps, function(tp) pred_sample_at(meta, s, tp), character(1)
    )
    if (anyNA(ids) || !all(ids %in% colnames(feature_mat))) {
      return(NULL)
    }
    if (spec$diff) {
      feature_mat[, ids[[2]]] - feature_mat[, ids[[1]]]
    } else {
      feature_mat[, ids[[1]]]
    }
  })
  names(rows) <- subjects
  rows <- rows[!vapply(rows, is.null, logical(1))]
  x <- do.call(rbind, rows)
  rownames(x) <- names(rows)
  x
}

# Subject-level outcome vector, named by subject. group = 1 for HR, 0 for LR;
# continuous outcomes come from phenotype.csv.
pred_outcome <- function(bundle, name) {
  if (name == "group") {
    grp <- bundle$meta |> distinct(subject, group)
    out <- as.integer(grp$group == "HR")
    names(out) <- grp$subject
    return(out)
  }
  ph <- bundle$pheno
  out <- ph[[name]]
  names(out) <- ph$subject
  out
}
