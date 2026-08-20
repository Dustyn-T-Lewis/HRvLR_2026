# The nine HRvLR contrasts and the pi-score threshold, shared by the non-imputed
# and imputed DEP runs so the two arms cannot drift.
#
# Each contrast is a linear combination of the six Group_Time cell means, which
# are estimable only for proteins observed in every cell they touch. The
# missingness filter keeps a protein detected in one cell alone (min_groups =
# 1), so 34 proteins reach the model with at least one empty cell. Their cell
# mean is inestimable, limma returns NA, and the protein is never tested. The
# true tested-N is 1877-1892 depending on the contrast, never 1900, and
# p.adjust sets the BH denominator from the non-NA count;
# b_reports/bh_denominators.csv records it per contrast.
# 03_Features/01_Proteins/a_script/03_untested_proteins.R names them.
HRVLR_CONTRASTS <- c(
  "Training_HR = HR_T2 - HR_T1",
  "Training_LR = LR_T2 - LR_T1",
  "Acute_HR = HR_T3 - HR_T2",
  "Acute_LR = LR_T3 - LR_T2",
  "Baseline_HRvLR = HR_T1 - LR_T1",
  "Trained_HRvLR = HR_T2 - LR_T2",
  "Acute_HRvLR = HR_T3 - LR_T3",
  "Training_Interaction = (HR_T2 - HR_T1) - (LR_T2 - LR_T1)",
  "Acute_Interaction = (HR_T3 - HR_T2) - (LR_T3 - LR_T2)"
)

CONTRAST_NAMES <- trimws(sub("=.*$", "", HRVLR_CONTRASTS))

# The continuous tree's two pooled contrasts: no arm term, averaged across all
# 16 subjects. Training_All and Acute_All were previously defined inline in
# 05_pooled_response.R; this is that same pair, promoted beside the nine so
# both trees read one contrasts file.
POOLED_CONTRASTS <- c(
  "Training = T2 - T1",
  "Acute = T3 - T2"
)

# proteoDA writes the fit wide, one column per statistic per contrast. The
# feature layer and the imputed arms read it long, so the pivot lives here
# instead of in each caller. Contrast names contain underscores and Acute_HR is
# a prefix of Acute_HRvLR, so the longest name has to be offered first.
dep_contrasts_long <- function(
  path = here::here(
    "03_Analysis", "categorical", "01_Proteins", "c_data",
    "03_combined_results.csv"
  )
) {
  ordered <- CONTRAST_NAMES[order(nchar(CONTRAST_NAMES), decreasing = TRUE)]
  pattern <- paste0("^(.*)_(", paste(ordered, collapse = "|"), ")$")
  readr::read_csv(path, show_col_types = FALSE) |>
    tidyr::pivot_longer(
      cols = tidyr::matches(pattern),
      names_pattern = pattern,
      names_to = c(".value", "contrast")
    ) |>
    dplyr::select(
      contrast, uniprot_id, gene, protein, description,
      logFC, CI.L, CI.R, average_intensity, t, B,
      P.Value, adj.P.Val, sig.PVal, sig.FDR, pi_score, sig_pi
    ) |>
    dplyr::arrange(match(contrast, CONTRAST_NAMES))
}

# The cell order the design matrix inherits. Pinned here because every arm has
# to place its columns identically for the contrast strings above to mean
# anything.
GROUP_LEVELS <- c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")

# The five contrasts that carry the HR-vs-LR question: three group differences
# and the two interactions. The four within-arm changes are the other half.
RESPONDER_CONTRASTS <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)

PI_THRESH <- 0.05

# Xiao et al. 2014 Eq. 2: Pi = p^|log2FC|, bounded in [0,1], lower = more
# significant. This is not the pi-value of Eq. 1 and carries no FDR control.
#
# sig_pi folds the gate and the direction into one column: +1 up, -1 down, 0
# otherwise. The 34 proteins limma could not test arrive with P.Value NA, so
# their pi_score is NA and both gated arms return NA; the final case_when
# branch is what lands them on 0L rather than NA. Untested and
# tested-but-unselected are deliberately indistinguishable here --
# 03_untested_proteins.R is what separates them.
pi_score <- function(p, logfc) p^abs(logfc)

add_pi_score <- function(res, pi_thresh = PI_THRESH) {
  res$pi_score <- pi_score(res$P.Value, res$logFC)
  res$sig_pi <- dplyr::case_when(
    res$pi_score < pi_thresh & res$logFC > 0 ~ 1L,
    res$pi_score < pi_thresh & res$logFC < 0 ~ -1L,
    TRUE ~ 0L
  )
  res
}
