# Phenotype mapping by association, not prediction. Each continuous training
# response (d_mcsa, d_1rm, ...) is one scalar per subject; we collapse the
# feature matrix to one value per subject per phase (baseline level, or the
# training/acute delta) and regress it on the trait with limma empirical-Bayes
# moderation. One observation per subject makes it leakage-free by construction,
# and moderation stabilises the per-feature variance at n<=16 (Smyth 2004;
# Ritchie 2015). Pathways use the same association on a singscore matrix, which
# is scored per sample and so cannot leak across folds (Foroutan 2018).

pacman::p_load(limma, dplyr, tibble, purrr, singscore)

# feature x sample matrix -> feature x subject matrix for one phase.
phase_subject_matrix <- function(mat, meta, phase,
                                 subject_col = "Subject_ID",
                                 time_col = "Timepoint") {
  tps <- switch(phase,
    baseline = "T1",
    training = c("T1", "T2"),
    acute = c("T2", "T3"),
    stop("phase must be baseline, training, or acute")
  )
  keep <- meta[[time_col]] %in% tps
  mat <- mat[, keep, drop = FALSE]
  subj <- meta[[subject_col]][keep]
  time <- meta[[time_col]][keep]
  if (phase == "baseline") {
    colnames(mat) <- subj
    return(mat)
  }
  early <- setNames(which(time == tps[1]), subj[time == tps[1]])
  late <- setNames(which(time == tps[2]), subj[time == tps[2]])
  both <- intersect(names(early), names(late))
  delta <- mat[, late[both], drop = FALSE] - mat[, early[both], drop = FALSE]
  colnames(delta) <- both
  delta
}

# limma association of every feature against one trait across subjects.
.fit_one_trait <- function(fmat, trait_values, trait, phase) {
  ok <- !is.na(trait_values)
  design <- model.matrix(~ trait_values[ok])
  fit <- limma::eBayes(limma::lmFit(fmat[, ok, drop = FALSE], design))
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tibble::tibble(
    feature = rownames(tt), trait = trait, phase = phase,
    beta = tt$logFC, t = tt$t, p = tt$P.Value, bh = tt$adj.P.Val
  )
}

# Associate every feature with every trait for one phase. `pheno` is one row per
# subject with a `subject` column; `traits` are its continuous-response columns.
associate_traits <- function(mat, meta, pheno, traits, phase,
                             subject_col = "Subject_ID") {
  fmat <- phase_subject_matrix(mat, meta, phase, subject_col)
  ph <- pheno[match(colnames(fmat), pheno$subject), , drop = FALSE]
  purrr::map_dfr(traits, \(tr) .fit_one_trait(fmat, ph[[tr]], tr, phase))
}

# Single-sample, rank-based pathway scores (pathways x samples). `gene_sets` is
# a named list of gene-symbol vectors; `mat` must be indexed by those symbols.
score_pathways <- function(mat, gene_sets, min_size = 5) {
  sizes <- vapply(gene_sets, \(g) sum(g %in% rownames(mat)), integer(1))
  gene_sets <- gene_sets[sizes >= min_size]
  ranked <- singscore::rankGenes(mat)
  singscore::multiScore(ranked, upSetColc = gene_sets)$Scores
}
