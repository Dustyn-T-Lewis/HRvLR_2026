# One estimator for every candidate label. Each label is a per-subject split
# into hi and lo, crossed with timepoint to give the same six cells the V1
# design used, so the two contrast strings below are literally the same for all
# seven labels and a logFC read off one sweep row means what it means on any
# other.
#
# The given HR/LR label is the identity case: its Baseline and
# Training_Interaction are two of V1's committed nine, and
# verify_v1_equivalence() refits them through this path to prove the
# generalisation did not move the numbers.

pacman::p_load(here, dplyr, tibble, readr, limma, openxlsx)
source(here("functions", "feature_levels.R"))

LABEL_CELLS <- c("hi_T1", "hi_T2", "hi_T3", "lo_T1", "lo_T2", "lo_T3")

LABEL_CONTRASTS <- c(
  "Baseline = hi_T1 - lo_T1",
  "Training_Interaction = (hi_T2 - hi_T1) - (lo_T2 - lo_T1)"
)

LABEL_CONTRAST_NAMES <- trimws(sub("=.*$", "", LABEL_CONTRASTS))

sample_metadata <- function() {
  dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
  m <- as.data.frame(dal$metadata)
  tibble(
    sample_id = m$Col_ID,
    subject = m$Subject_ID,
    timepoint = m$Timepoint
  )
}

# `label` is named by subject and holds "hi"/"lo". Subjects it does not name
# are dropped, which is how a label built on d_1rm_ext loses HR_S22 to its one
# NA without the caller having to special-case it.
label_cells <- function(label, meta = sample_metadata()) {
  meta |>
    filter(.data$subject %in% names(label)) |>
    mutate(cell = factor(
      paste0(unname(label[.data$subject]), "_", .data$timepoint),
      levels = LABEL_CELLS
    ))
}

label_design <- function(mat, label, adjust = NULL,
                         meta = sample_metadata()) {
  meta <- label_cells(label, meta)
  meta <- meta[match(colnames(mat), meta$sample_id), ]
  keep <- !is.na(meta$sample_id)
  mat <- mat[, keep, drop = FALSE]
  meta <- meta[keep, ]
  empty <- setdiff(LABEL_CELLS, as.character(meta$cell))
  if (length(empty)) {
    stop("label leaves no samples in cell(s): ", paste(empty, collapse = ", "))
  }
  design <- stats::model.matrix(~ 0 + cell, meta)
  colnames(design) <- levels(meta$cell)
  if (!is.null(adjust)) {
    design <- cbind(design, adjust = unname(adjust[colnames(mat)]))
  }
  cm <- limma::makeContrasts(contrasts = LABEL_CONTRASTS, levels = design)
  colnames(cm) <- LABEL_CONTRAST_NAMES
  list(
    mat = mat, design = design, contrasts = cm, block = meta$subject,
    correlation = limma::duplicateCorrelation(
      mat, design,
      block = meta$subject
    )$consensus
  )
}

# BH within each contrast, never pooled across contrasts or across labels. The
# seven labels overlap heavily by construction -- five are splits of correlated
# outcomes on the same 16 subjects -- so a q spanning them would claim an
# independence the design does not have. The sweep reports its own cell count
# instead.
fit_label_contrasts <- function(mat, label, adjust = NULL, robust = FALSE) {
  parts <- label_design(mat, label, adjust = adjust)
  fit <- limma::lmFit(parts$mat, parts$design,
    block = parts$block, correlation = parts$correlation
  )
  fit2 <- limma::eBayes(
    limma::contrasts.fit(fit, parts$contrasts),
    robust = robust
  )
  res <- bind_rows(lapply(LABEL_CONTRAST_NAMES, function(ct) {
    limma::topTable(fit2,
      coef = ct, number = Inf, adjust.method = "BH", sort.by = "none"
    ) |>
      tibble::rownames_to_column("feature") |>
      transmute(
        contrast = ct, feature = .data$feature, logFC = .data$logFC,
        t = .data$t, p = .data$P.Value, bh = .data$adj.P.Val
      )
  }))
  attr(res, "within_cor") <- parts$correlation
  attr(res, "n_samples") <- ncol(parts$mat)
  res
}

# The generalisation check: Baseline and Training_Interaction under the given
# HR/LR label are V1's Baseline_HRvLR and Training_Interaction. robust = TRUE
# matches proteoDA's own eBayes call, which is what the V1 workbook was
# written from.
V1_RESULTS <- file.path(
  dirname(here()), "A_HRvLR_2026", "03_Analysis", "categorical",
  "01_Proteins", "c_data", "05_results.xlsx"
)

verify_v1_equivalence <- function(book = V1_RESULTS, tol = 1e-6) {
  ph <- read_csv(here("00_input", "c_data", "phenotype.csv"),
    show_col_types = FALSE
  )
  given <- setNames(ifelse(ph$group_arm == "HR", "hi", "lo"), ph$subject)
  fitted <- fit_label_contrasts(protein_matrix(), given, robust = TRUE)
  v1_name <- c(
    Baseline = "Baseline_HRvLR",
    Training_Interaction = "Training_Interaction"
  )
  bind_rows(lapply(LABEL_CONTRAST_NAMES, function(ct) {
    ref <- openxlsx::read.xlsx(book, v1_name[[ct]])
    j <- inner_join(
      filter(fitted, .data$contrast == ct),
      transmute(ref,
        feature = .data$uniprot_id, ref_lfc = .data$logFC,
        ref_p = .data$P.Value, ref_bh = .data$adj.P.Val
      ),
      by = "feature"
    )
    worst <- max(
      abs(j$logFC - j$ref_lfc), abs(j$p - j$ref_p), abs(j$bh - j$ref_bh),
      na.rm = TRUE
    )
    tibble(
      contrast = ct, v1_contrast = v1_name[[ct]], n = nrow(j),
      max_abs_diff = worst, equivalent = worst < tol
    )
  }))
}
