# What the HR/LR label is made of, and what it does and does not separate.
#
# The label is downstream of comp_hypertrophy, which arrives from the
# collaborator's spreadsheet without a stated formula. This stage recovers what
# it behaves like, records whether each outcome changed over training at all,
# and reports the separation the label achieves on each. It draws no inference
# about the V1 proteomic null; that belongs to the discussion.

pacman::p_load(
  here, dplyr, tidyr, purrr, tibble, readr, broom, mclust, openxlsx
)

OUT_DIR <- here("03_Features", "01_Responsiveness", "c_data")

TRAITS <- c(
  "d_fcsa_I", "d_fcsa_II", "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
)

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
label_separation <- map_dfr(c("comp_hypertrophy", TRAITS), function(v) {
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
composite_structure <- map_dfr(TRAITS, function(v) {
  fit <- stats::lm(pheno$comp_hypertrophy ~ pheno[[v]])
  tibble(
    trait = v,
    r = cor(pheno$comp_hypertrophy, pheno[[v]], use = "complete.obs"),
    r2_alone = summary(fit)$r.squared
  )
})

joint_r2 <- summary(stats::lm(
  comp_hypertrophy ~ d_fcsa_I + d_fcsa_II + d_mcsa +
    d_1rm_legpress + d_1rm_ext,
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
set.seed(42)
mod <- mclust::Mclust(pheno$comp_hypertrophy, G = 1:4, verbose = FALSE)

# BIC alone is thin at n = 16, and a two-component 1D mixture will always find
# the largest gap whether or not it means anything. The bootstrap LRT is the
# test that carries a p-value: it simulates from the fitted G-component model
# to build the null distribution of the likelihood ratio against G + 1.
lrt <- mclust::mclustBootstrapLRT(
  pheno$comp_hypertrophy,
  modelName = mod$modelName, maxG = 2, nboot = 999
)

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

# Two fibre-CSA measurement methods on the same samples. Reported because they
# disagree in a way that cannot be a scale difference, not because anything
# downstream uses them.
meta <- read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE)
method_check <- map_dfr(
  list(
    c("fCSA_Mixed_Pre", "MyoVision_fCSA_mixed_Pre"),
    c("fCSA_Type_I_Pre", "MyoVision_fCSA_Type_I__Pre")
  ),
  function(pair) {
    d <- meta |>
      filter(.data$Timepoint %in% c("T1", "T2")) |>
      select(manual = all_of(pair[1]), myovision = all_of(pair[2])) |>
      tidyr::drop_na()
    tibble(
      quantity = pair[1], n = nrow(d),
      r = cor(d$manual, d$myovision),
      mean_difference = mean(d$myovision - d$manual),
      sd_difference = sd(d$myovision - d$manual)
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
    method_check = method_check
  ),
  file.path(OUT_DIR, "01_label_audit.xlsx")
)

message(
  "composite joint r2 on the five outcomes: ", round(joint_r2, 3),
  "; label is the exact median cut: ",
  split_check$top_half_all_hr && split_check$bottom_half_all_lr
)
