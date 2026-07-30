# The nine HRvLR contrasts and the pi-score threshold, shared by the non-imputed
# and imputed DEP runs so the two arms cannot drift.
#
# Each contrast is a linear combination of the six Group_Time cell means, which are estimable
# only for proteins observed in every cell they touch. The missingness filter keeps a protein
# detected in one cell alone (min_groups = 1), so 34 proteins reach the model with at least one
# empty cell. Their cell mean is inestimable, limma returns NA, and the protein is never tested.
# The true tested-N is 1882-1897 depending on the contrast, never 1905, and p.adjust sets the BH
# denominator from the non-NA count; b_reports/bh_denominators.csv records it
# per contrast. 03_DEP/a_non_imputed/a_script/03_untested_proteins.R names them.
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

PI_THRESH <- 0.05

# Xiao et al. 2014 Eq. 2: Pi = p^|log2FC|, bounded in [0,1], lower = more significant.
# This is not the pi-value of Eq. 1 and carries no FDR control.
#
# sig_pi folds the gate and the direction into one column: +1 up, -1 down, 0 otherwise.
# The 34 proteins limma could not test arrive with P.Value NA, so their pi_score is NA and
# both gated arms return NA; the final case_when branch is what lands them on 0L rather
# than NA. Untested and tested-but-unselected are deliberately indistinguishable here --
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
