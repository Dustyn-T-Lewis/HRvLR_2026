# Shared inputs and feature/outcome builders for the prediction suite.
# Loads the imputed proteome, collapses to genes, scores singscore pathways
# once on the full cohort, reshapes the F04_association/WGCNA module eigengenes, and exposes
# per-phase subject-level feature matrices plus the phenotype outcomes.
#
# singscore is single-sample and rank-based: each sample's pathway score
# depends only on that sample's own within-sample gene ranks, never on the
# rest of the cohort. Scoring the full matrix once is therefore leakage-free,
# which is why the suite uses singscore rather than a cohort-relative method
# such as GSVA. Every fold-specific transform (z-scoring, model tuning) still
# happens train-only inside the harness.

pacman::p_load(here, dplyr, tidyr, readr, tibble, limma, singscore)
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_singscore.R"))

pred_paths <- function() {
  list(
    dalist = here(
      "02_Normalization", "imputation", "c_data",
      "DAList_imputed_missforest.rds"
    ),
    pheno = here("00_input", "c_data", "phenotype.csv"),
    eigen = here("04_Figures", "F04_association", "WGCNA", "c_data", "wgcna_eigengene.csv"),
    cache = here(
      "04_Figures", "F06_prediction", "prediction_shared", "c_data",
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

# Load everything once and return the shared bundle.
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
    singscore  = pred_singscore(expr, paths$cache),
    proteins   = as.matrix(da$data),
    eigengenes = pred_eigengene_matrix(paths$eigen, colnames(da$data))
  )

  pheno <- read_csv(paths$pheno, show_col_types = FALSE)

  list(meta = meta, feature_sets = feature_sets, pheno = pheno)
}

# Sample id for a given subject/timepoint, or NA when that sample is absent.
pred_sample_at <- function(meta, subj, tp) {
  hit <- meta$sample[meta$subject == subj & meta$timepoint == tp]
  if (length(hit) == 1) hit else NA_character_
}

# Subject-level feature matrix for one phase: baseline is the T1 sample,
# training is the T2 - T1 change, acute is the T3 - T2 change. Subjects
# missing a required sample are dropped.
pred_phase_matrix <- function(feature_mat, meta, phase) {
  tp_needed <- switch(phase,
    baseline = "T1",
    training = c("T1", "T2"),
    acute    = c("T2", "T3"),
    stop("unknown phase: ", phase)
  )
  subjects <- levels(factor(meta$subject))

  rows <- lapply(subjects, function(s) {
    ids <- vapply(
      tp_needed, function(tp) pred_sample_at(meta, s, tp),
      character(1)
    )
    if (anyNA(ids) || !all(ids %in% colnames(feature_mat))) {
      return(NULL)
    }
    if (phase == "baseline") {
      feature_mat[, ids[[1]]]
    } else {
      feature_mat[, ids[[2]]] - feature_mat[, ids[[1]]]
    }
  })
  names(rows) <- subjects
  rows <- rows[!vapply(rows, is.null, logical(1))]

  x <- do.call(rbind, rows)
  rownames(x) <- names(rows)
  x
}

# Subject-level outcome vector, named by subject. For the class arm the
# outcome is 1 = HR, 0 = LR; continuous outcomes come from phenotype.csv.
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
