# What the HR/LR label is made of, and what it does and does not separate.
#
# The label is downstream of comp_hypertrophy, which arrives from the
# collaborator's spreadsheet without a stated formula. This stage recovers what
# it behaves like, records whether each outcome changed over training at all,
# and reports the separation the label achieves on each. It draws no inference
# about the V1 proteomic null; that belongs to the discussion.

pacman::p_load(
  here, dplyr, tidyr, purrr, tibble, readr, broom, withr, openxlsx
)

# mclust is attached per call rather than for the session. Mclust() evaluates
# its own matched call in the caller's frame, so it needs the package on the
# search path and cannot be used purely qualified; leaving it attached masks
# purrr::map with mclust::map for every script sourced afterwards.

OUT_DIR <- here("03_Features", "01_Responsiveness", "c_data")

# Change scores. volume_load is excluded because it is a total, not a change:
# asking whether it differs from zero is asking whether the cohort lifted
# anything. It still gets a candidate label and still appears in the separation
# table, where the question is meaningful.
TRAITS <- c(
  "d_fcsa_I", "d_fcsa_II", "d_fcsa_mixed", "d_nfibre_mixed", "d_nfibre_I",
  "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
)

SPLIT_TRAITS <- c(TRAITS, "volume_load")

pheno <- read_csv(here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
)

# Did the cohort change at all? A trait whose change CI spans zero cannot
# support a responder ordering no matter how the split is drawn.
change_summary <- map_dfr(c("comp_hypertrophy", TRAITS), function(v) {
  broom::tidy(stats::t.test(pheno[[v]])) |>
    transmute(
      trait = v, n = sum(!is.na(pheno[[v]])), mean = estimate,
      ci_lo = conf.low, ci_hi = conf.high, t = statistic, p = p.value
    )
}) |>
  mutate(
    sd_change = map_dbl(trait, ~ sd(pheno[[.x]], na.rm = TRUE)),
    across(c(mean, ci_lo, ci_hi), ~ .x / sd_change, .names = "{.col}_d")
  )

# What the label separates. comp_hypertrophy is included as the identity case:
# the label is its median split, so its separation is arithmetic, not evidence.
label_separation <- map_dfr(c("comp_hypertrophy", SPLIT_TRAITS), function(v) {
  broom::tidy(stats::t.test(pheno[[v]] ~ pheno$group_arm)) |>
    transmute(
      trait = v, hr = estimate1, lr = estimate2,
      difference = estimate1 - estimate2,
      ci_lo = conf.low, ci_hi = conf.high,
      t = statistic, p = p.value
    )
}) |>
  mutate(
    sd_change = map_dbl(trait, ~ sd(pheno[[.x]], na.rm = TRUE)),
    across(c(difference, ci_lo, ci_hi), ~ .x / sd_change, .names = "{.col}_d")
  )

# How much of comp_hypertrophy each outcome accounts for. The single-predictor
# r2 column is what decides the internal/external flag the sweep carries.
composite_structure <- map_dfr(SPLIT_TRAITS, function(v) {
  fit <- stats::lm(pheno$comp_hypertrophy ~ pheno[[v]])
  tibble(
    trait = v,
    r = cor(pheno$comp_hypertrophy, pheno[[v]], use = "complete.obs"),
    r2_alone = summary(fit)$r.squared
  )
})

joint_r2 <- summary(stats::lm(
  stats::reformulate(TRAITS, response = "comp_hypertrophy"),
  data = pheno
))$r.squared

# Is the label exactly the median cut of the composite?
ranked <- pheno |> arrange(desc(comp_hypertrophy))
split_check <- tibble(
  n = nrow(ranked),
  top_half_all_hr = all(ranked$group_arm[1:8] == "HR"),
  bottom_half_all_lr = all(ranked$group_arm[9:16] == "LR"),
  gap_at_cut = ranked$comp_hypertrophy[8] - ranked$comp_hypertrophy[9],
  sd_composite = sd(pheno$comp_hypertrophy)
)

# The cut sits in a wide gap, so ask whether the composite is actually bimodal
# rather than a median split through a continuum. mclust selects the component
# count by BIC over G = 1:4 and can return G = 1, which is the answer that says
# there is no grouping to find. This is the phenotype-side counterpart of the
# question stage 03 asks of the proteome.
# BIC alone is thin at n = 16, and a two-component 1D mixture will always find
# the largest gap whether or not it means anything. The bootstrap LRT is the
# test that carries a p-value: it simulates from the fitted G-component model
# to build the null distribution of the likelihood ratio against G + 1.
composite_mixture <- function(x) {
  withr::local_package("mclust")
  fit <- mclust::Mclust(x, G = 1:4, verbose = FALSE)
  list(
    fit = fit,
    lrt = mclust::mclustBootstrapLRT(
      x,
      modelName = fit$modelName, maxG = 2, nboot = 999
    )
  )
}

set.seed(42)
fitted_mixture <- composite_mixture(pheno$comp_hypertrophy)
mod <- fitted_mixture$fit
lrt <- fitted_mixture$lrt

composite_modality <- tibble(
  best_g = mod$G,
  model = mod$modelName,
  bic_g1 = mod$BIC[1, 1],
  bic_best = max(mod$BIC, na.rm = TRUE),
  bic_gain_over_g1 = max(mod$BIC, na.rm = TRUE) - mod$BIC[1, 1],
  lrt_g1_vs_g2 = lrt$obs[1],
  lrt_p = lrt$p.value[1],
  agreement_with_label = mclust::adjustedRandIndex(
    mod$classification, pheno$group_arm
  )
)

# The MyoVision columns are fibre counts, not areas, despite the "fCSA" in
# their meta names; the source workbook calls them "Number of fCSA - Mixed
# (MyoVision)". Recorded here because the naming invites reading them as a
# second area measurement that disagrees with the first, and they do not
# disagree: a count runs an order of magnitude below an area and moves against
# it, since larger fibres pack fewer into the imaged field.
meta <- read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE)
fibre_count_check <- map_dfr(
  list(
    c("fCSA_Mixed_Pre", "MyoVision_fCSA_mixed_Pre"),
    c("fCSA_Type_I_Pre", "MyoVision_fCSA_Type_I__Pre")
  ),
  function(pair) {
    d <- meta |>
      filter(.data$Timepoint %in% c("T1", "T2")) |>
      select(area = all_of(pair[1]), count = all_of(pair[2])) |>
      tidyr::drop_na()
    tibble(
      quantity = pair[1], n = nrow(d),
      area_median = stats::median(d$area),
      count_median = stats::median(d$count),
      r_area_count = cor(d$area, d$count),
      r_area_inverse_count = cor(d$area, 1 / d$count)
    )
  }
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.xlsx(
  list(
    change_summary = change_summary,
    label_separation = label_separation,
    composite_structure = composite_structure,
    split_check = split_check,
    composite_modality = composite_modality,
    fibre_count_check = fibre_count_check
  ),
  file.path(OUT_DIR, "01_label_audit.xlsx")
)

message(
  "composite joint r2 on the eight change scores: ", round(joint_r2, 3),
  "; label is the exact median cut: ",
  split_check$top_half_all_hr && split_check$bottom_half_all_lr
)
