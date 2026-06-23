# Shared input loader for 05_Clustering. All engine scripts source this file
# and call load_clustering_inputs() to get a consistent, pi-gated input set.
# pi-gated = pi_score < 0.05 in ANY contrast (Xiao π-value, P^|logFC|).
# Abundance is from the canonical missforest-imputed matrix (zero NAs).
# Gap matrix columns: HR−LR logFC at T1 (Baseline), T2 (Trained), T3 (Acute).

pacman::p_load(here, dplyr)

load_clustering_inputs <- function() {
  primary <- readRDS(here("03_DEP", "c_data", "01_limma_DAList.rds"))
  imputed <- readRDS(here(
    "02_Normalization", "imputation", "c_data",
    "DAList_imputed_missforest.rds"
  ))
  stopifnot(identical(imputed$metadata$Col_ID, colnames(imputed$data)))

  # pi-gate: proteins with pi_score < 0.05 in any contrast
  pi_mat <- vapply(
    primary$results, `[[`, numeric(nrow(primary$data)), "pi_score"
  )
  rownames(pi_mat) <- rownames(primary$data)
  pi_keep <- apply(pi_mat, 1, function(x) any(x < 0.05, na.rm = TRUE))
  pi_set <- rownames(pi_mat)[pi_keep]
  n_before <- length(pi_set)
  pi_set <- intersect(pi_set, rownames(imputed$data))
  message(
    n_before - length(pi_set),
    " pi-gated proteins dropped (absent from imputed matrix)"
  )

  abund <- imputed$data[pi_set, , drop = FALSE]

  md <- imputed$metadata
  meta <- data.frame(
    sample_id = md$Col_ID,
    subject = md$Subject_ID,
    group_arm = md$Group,
    timepoint = md$Timepoint,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # gap: HR−LR logFC at each timepoint, rows restricted to pi_set
  hrvlr_contrasts <- c("Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR")
  gap <- vapply(
    primary$results[hrvlr_contrasts],
    function(r) r$logFC, numeric(nrow(primary$data))
  )
  rownames(gap) <- rownames(primary$data)
  gap <- gap[pi_set, , drop = FALSE]
  colnames(gap) <- c("T1", "T2", "T3")

  stopifnot(
    ncol(abund) == 45,
    ncol(gap) == 3,
    nrow(abund) == length(pi_set),
    !anyNA(abund),
    identical(meta$sample_id, colnames(abund))
  )

  list(abund = abund, meta = meta, gap = gap, pi_set = pi_set)
}
