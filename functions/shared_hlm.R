# Mixed-model association fitter: builds the standard response ~ fixed +
# (1 | group) formula and returns the fitted lmerTest model, so every caller
# shares one model spec and keeps its own fixef/Satterthwaite extraction.
fit_hlm <- function(data, response, fixed = "Group * Timepoint",
                    group = "subject", reml = TRUE) {
  pacman::p_load(lmerTest)
  form <- stats::as.formula(
    sprintf("%s ~ %s + (1 | %s)", response, fixed, group)
  )
  lmerTest::lmer(form, data = data, REML = reml)
}
