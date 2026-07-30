# Did both arms get equally bloodier at T3?
#
# T3 biopsies carry roughly twice the haemoglobin signal of T1 and T2. F04 says
# the confound cancels in the interaction contrasts, which holds only if the
# rise is the same size in HR and LR. Nothing had tested that. Filtering removes
# blood proteins, not blood, so a residual arm difference in carryover would sit
# underneath every acute contrast whatever the protein list looks like.
#
# The arm-by-T3 term is the answer, and it is reported whichever way it lands.

pacman::p_load(here, dplyr, tibble)
source(here("functions", "shared_hlm.R"))

# The three metadata rows the pipeline excluded (the partial-roster biopsies)
# carry a blood index but no analysed sample, so they are dropped here too.
# Fitting all 48 would model biopsies that no contrast ever saw.
blood_index_data <- function() {
  intermediates <- readRDS(
    here("01_Filtering", "c_data", "filtering_intermediates.rds")
  )
  analysed <- colnames(
    readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))$data
  )
  intermediates$blood_index |>
    filter(.data$Col_ID %in% analysed) |>
    transmute(
      Col_ID = .data$Col_ID,
      subject = .data$Subject_ID,
      Group = factor(.data$Group, levels = c("LR", "HR")),
      Timepoint = factor(.data$Timepoint, levels = c("T1", "T2", "T3")),
      blood_index = .data$blood_index
    )
}

blood_index_arm_test <- function(data = blood_index_data()) {
  co <- summary(fit_hlm(data, "blood_index"))$coefficients
  tibble(
    term = rownames(co),
    estimate = unname(co[, "Estimate"]),
    se = unname(co[, "Std. Error"]),
    df = unname(co[, "df"]),
    t = unname(co[, "t value"]),
    p = unname(co[, "Pr(>|t|)"])
  )
}

# The one line F04 prints. Satterthwaite df are fractional, so they are rounded
# for display only.
blood_arm_t3_label <- function(terms = blood_index_arm_test()) {
  row <- terms[terms$term == "GroupHR:TimepointT3", ]
  sprintf(
    "arm x T3 on the blood index: b = %+.2f (SE %.2f, df %.0f, p = %.3f)",
    row$estimate, row$se, row$df, row$p
  )
}
