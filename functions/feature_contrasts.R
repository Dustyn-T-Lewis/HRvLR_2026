# One estimator for all three feature levels. Stage 03 fits the nine contrasts
# on proteins through proteoDA; modules and pathways are fitted here on the same
# means model with the same duplicateCorrelation blocking, so a logFC read off a
# module row and a logFC read off a protein row mean the same thing.
# verify_protein_equivalence() is the proof: it refits the protein matrix
# through this path and checks it against the committed stage 03 workbook.

pacman::p_load(here, dplyr, tibble, readr, limma, openxlsx)
source(here("03_DEP", "contrasts.R"))

GROUP_LEVELS <- c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")

feature_metadata <- function() {
  dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
  m <- as.data.frame(dal$metadata)
  tibble(
    sample_id = m$Col_ID,
    group = factor(m$Group_Time, levels = GROUP_LEVELS),
    subject = m$Subject_ID
  )
}

# BH is applied within each contrast, matching extract_DA_results()'s per-coef
# topTable call. Pooling the nine would imply an independence that the shared
# subjects and overlapping contrasts do not have.
fit_feature_contrasts <- function(mat, meta = feature_metadata()) {
  meta <- meta[match(colnames(mat), meta$sample_id), ]
  if (anyNA(meta$sample_id)) {
    stop("feature matrix has samples absent from the normalisation metadata")
  }
  design <- stats::model.matrix(~ 0 + group, meta)
  colnames(design) <- levels(meta$group)
  within_cor <- limma::duplicateCorrelation(mat, design, block = meta$subject)
  fit <- limma::lmFit(mat, design,
    block = meta$subject, correlation = within_cor$consensus
  )
  cm <- limma::makeContrasts(contrasts = HRVLR_CONTRASTS, levels = design)
  colnames(cm) <- trimws(sub("=.*$", "", HRVLR_CONTRASTS))
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
  res <- bind_rows(lapply(colnames(cm), function(ct) {
    limma::topTable(fit2,
      coef = ct, number = Inf, adjust.method = "BH", sort.by = "none"
    ) |>
      tibble::rownames_to_column("feature") |>
      transmute(
        contrast = ct, feature = .data$feature, logFC = .data$logFC,
        t = .data$t, p = .data$P.Value, bh = .data$adj.P.Val
      )
  }))
  attr(res, "within_cor") <- within_cor$consensus
  res
}

protein_matrix <- function() {
  df <- read_csv(here("02_Normalization", "c_data", "normalized.csv"),
    show_col_types = FALSE
  )
  ann_cols <- c("uniprot_id", "protein", "gene", "description")
  mat <- as.matrix(df[, setdiff(names(df), ann_cols)])
  rownames(mat) <- df$uniprot_id
  mat
}

verify_protein_equivalence <- function(tol = 1e-6) {
  fitted <- fit_feature_contrasts(protein_matrix())
  book <- here("03_DEP", "a_non_imputed", "c_data", "05_results.xlsx")
  bind_rows(lapply(unique(fitted$contrast), function(ct) {
    ref <- openxlsx::read.xlsx(book, ct)
    j <- inner_join(
      filter(fitted, .data$contrast == ct),
      transmute(ref,
        feature = .data$uniprot_id, ref_lfc = .data$logFC,
        ref_p = .data$P.Value, ref_bh = .data$adj.P.Val
      ),
      by = "feature"
    )
    tibble(
      contrast = ct, n = nrow(j),
      max_d_logFC = max(abs(j$logFC - j$ref_lfc), na.rm = TRUE),
      max_d_p = max(abs(j$p - j$ref_p), na.rm = TRUE),
      max_d_bh = max(abs(j$bh - j$ref_bh), na.rm = TRUE),
      equivalent = max(
        abs(j$logFC - j$ref_lfc), abs(j$p - j$ref_p), abs(j$bh - j$ref_bh),
        na.rm = TRUE
      ) < tol
    )
  }))
}
