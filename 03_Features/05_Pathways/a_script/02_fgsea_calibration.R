# Is fgsea's padj trustworthy on this proteome?
#
# 01_run_pathways.R returned 1398 set-contrast hits at BH < 0.05, down to
# padj = 2e-19, under labels whose protein-level fits produced nothing that
# survived a permutation check. It also returned more hits under 1RM leg
# extension, a split with essentially no agreement with the responder label,
# than under the responder label itself. Meanwhile singscore, scored per sample
# and fitted through the estimator the proteins used, returned zero.
#
# fgsea's preranked null permutes gene labels, which treats proteins as
# exchangeable. Members of a pathway are co-regulated and, in a normalised MS
# matrix, share technical structure too, so that null is anticonservative here.
# This repo has already recorded the same failure mode for STRING's PPI
# enrichment p on this proteome.
#
# So the question is not which pathways came out but whether any count means
# anything. Shuffle the label across subjects, refit, rerun fgsea, and count.
# If a random split yields a comparable number of significant sets, the padj
# column carries no information and no pathway hit from stage 05 can stand.

pacman::p_load(here, dplyr, purrr, tibble, readr, limma, fgsea, openxlsx)

source(here("functions", "label_contrasts.R"))
source(here("functions", "shared_pathway_utils.R"))

OUT_DIR <- here("03_Features", "05_Pathways", "c_data")
BH_ALPHA <- 0.05
N_PERM <- 100
CONTRAST <- "Baseline"
LABEL <- "given"

labels_long <- read_csv(
  here("03_Features", "01_Responsiveness", "c_data", "02_candidate_labels.csv"),
  show_col_types = FALSE
)
annotation <- read_csv(
  here("02_Normalization", "c_data", "normalized.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(feature = uniprot_id, gene)

pathways <- build_pathway_collection()
mat <- protein_matrix()

lab <- filter(labels_long, .data$label == LABEL)
vec <- setNames(lab$level, lab$subject)
parts <- label_design(mat, vec)
samples <- colnames(parts$mat)

# One label assignment through to a hit count. The correlation and block are
# held fixed for the same reason 02_confirm_hits.R holds them: permuting who is
# hi does not disturb the repeated-measures structure they describe.
count_hits <- function(label_vec) {
  meta <- label_cells(label_vec)
  meta <- meta[match(samples, meta$sample_id), ]
  design <- stats::model.matrix(~ 0 + cell, meta)
  colnames(design) <- levels(meta$cell)
  cm <- limma::makeContrasts(contrasts = LABEL_CONTRASTS, levels = design)
  colnames(cm) <- LABEL_CONTRAST_NAMES
  fit <- limma::eBayes(limma::contrasts.fit(
    limma::lmFit(mat[, samples, drop = FALSE], design,
      block = meta$subject, correlation = parts$correlation
    ), cm
  ))
  tt <- limma::topTable(fit,
    coef = CONTRAST, number = Inf, sort.by = "none"
  ) |>
    tibble::rownames_to_column("feature")
  ranks <- tt |>
    left_join(annotation, by = "feature") |>
    filter(!is.na(.data$gene), .data$gene != "", !is.na(.data$t)) |>
    slice_max(abs(.data$t), n = 1, by = "gene", with_ties = FALSE) |>
    (\(d) setNames(d$t, d$gene))()
  res <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = pathways, stats = ranks,
    minSize = 15, maxSize = 500, nPermSimple = 1000, eps = 0
  ))
  sum(res$padj < BH_ALPHA, na.rm = TRUE)
}

observed <- count_hits(vec)

set.seed(42)
null_counts <- vapply(seq_len(N_PERM), function(i) {
  count_hits(setNames(sample(unname(vec)), names(vec)))
}, numeric(1))

calibration <- tibble(
  label = LABEL, contrast = CONTRAST, n_perm = N_PERM,
  observed_hits = observed,
  null_hits_median = stats::median(null_counts),
  null_hits_q05 = stats::quantile(null_counts, 0.05),
  null_hits_q95 = stats::quantile(null_counts, 0.95),
  null_any_hit = mean(null_counts > 0),
  p_empirical = (sum(null_counts >= observed) + 1) / (N_PERM + 1)
)

write_csv(
  tibble(perm = seq_len(N_PERM), n_hits = null_counts),
  file.path(OUT_DIR, "02_fgsea_null_counts.csv")
)
write.xlsx(
  list(calibration = calibration),
  file.path(OUT_DIR, "02_fgsea_calibration.xlsx")
)

print(as.data.frame(calibration), digits = 3)
message(
  "\nobserved ", observed, " significant sets; a random split of the same ",
  "subjects gives a median of ", stats::median(null_counts),
  " (empirical p = ", signif(calibration$p_empirical, 2), ")"
)
