# The nine HRvLR contrasts and the pi-score threshold, shared by the non-imputed
# and imputed DEP runs so the two arms cannot drift.
#
# Each contrast is a linear combination of the six Group_Time cell means, which are estimable
# only for proteins observed in every cell they touch. The missingness filter keeps a protein
# detected in one cell alone (min_groups = 1), so 34 proteins reach the model with at least one
# empty cell. Their cell mean is inestimable, limma returns NA, and the protein is never tested.
# The true tested-N is 1897-1912 depending on the contrast, never 1920, and p.adjust sets the BH
# denominator from the non-NA count. 03_DEP/a_non_imputed/a_script/03_untested_proteins.R names
# them.
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
