# Phenotype stats for F01: the measured outcomes, the mixed models, and the
# divergence estimate. No plotting - each panel owns its own geometry.
pacman::p_load(
  here, dplyr, tidyr, readr, tibble, purrr, forcats,
  effectsize, lmerTest, emmeans
)
source(here::here("functions", "shared_hlm.R"))

# Repeated-measures outcomes (T1 -> T2). The MyoVision_fCSA_* columns are left
# out: at 200-560 they sit 10-30x below plausible fibre CSA and correlate
# negatively with manual tracing, the signature of a fibre-count metric
# mislabeled as area (MyoVision reports CSA in thousands of um^2; Wen et al.
# 2018, doi:10.1152/japplphysiol.00762.2017).
MEASURES <- tibble::tribble(
  ~col, ~label, ~domain,
  "mCSA_Pre", "mCSA", "muscle",
  "fCSA_Mixed_Pre", "fCSA Mixed", "fibre",
  "fCSA_Type_I_Pre", "fCSA Type I", "fibre",
  "fCSA_Type_II_Pre", "fCSA Type II", "fibre",
  "X1RM_Leg_Pre", "1RM Leg Press", "strength",
  "X1RM._Ext_Pre", "1RM Extension", "strength"
)

f01_meta <- function() {
  read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE) |>
    mutate(
      subject = sub("_T[123]$", "", Col_ID),
      Group = factor(Group, levels = c("HR", "LR")),
      Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
    )
}

f01_composite_scores <- function(meta) {
  meta |>
    filter(Timepoint == "T2", !is.na(COMP.HYPERTROPHY)) |>
    transmute(subject, Group, value = readr::parse_number(COMP.HYPERTROPHY))
}

prepost_long <- function(meta, col) {
  meta |>
    filter(Timepoint %in% c("T1", "T2"), !is.na(.data[[col]])) |>
    transmute(
      subject, Group,
      Timepoint = factor(Timepoint, levels = c("T1", "T2")),
      value = .data[[col]]
    )
}

# The raw mixed model for one outcome, reused by the report tables and the
# diagnostics. Fitted on the response scale by REML with Satterthwaite df.
fit_lmm <- function(meta, col) {
  fit_hlm(prepost_long(meta, col), response = "value")
}

# The divergence estimate: each subject's T2-T1 change, then the HR-minus-LR
# standardized difference (Hedges g, >0 = HR gains more) with its 95% CI. At two
# timepoints this equals the Group x Timepoint mixed-model interaction (Twisk et
# al. 2018, doi:10.1016/j.conctc.2018.03.008); the hlm supplement shows the mixed
# model agrees. Fibre, muscle, and strength share the SD axis.
change_advantage <- function(meta, col) {
  w <- prepost_long(meta, col) |>
    pivot_wider(names_from = Timepoint, values_from = value) |>
    filter(!is.na(T1), !is.na(T2)) |>
    mutate(delta = T2 - T1)
  g <- effectsize::hedges_g(delta ~ Group, data = w)
  tibble(
    estimate = g$Hedges_g, ci_lo = g$CI_low, ci_hi = g$CI_high,
    p = t.test(delta ~ Group, data = w)$p.value
  )
}

# Every outcome's divergence, Holm-adjusted across the six-outcome family.
change_advantage_table <- function(meta) {
  purrr::map_dfr(seq_len(nrow(MEASURES)), function(i) {
    change_advantage(meta, MEASURES$col[i]) |>
      mutate(measure = MEASURES$label[i], domain = MEASURES$domain[i], .before = 1)
  }) |>
    mutate(p_holm = p.adjust(p, method = "holm"))
}

# The six raw models, named by outcome, for the easystats report and diagnostics.
lmm_models <- function(meta) {
  setNames(lapply(MEASURES$col, function(col) fit_lmm(meta, col)), MEASURES$label)
}

# Leave-one-subject influence on the divergence coefficient. A balanced 2x2 design
# gives every observation identical leverage, so check_model's leverage view is
# uninformative; this refit-based shift is the meaningful influence check.
lmm_loso_influence <- function(meta) {
  bind_rows(lapply(seq_len(nrow(MEASURES)), function(i) {
    d <- prepost_long(meta, MEASURES$col[i])
    d$z <- as.numeric(scale(d$value))
    beta <- function(x) lme4::fixef(fit_hlm(x, response = "z"))[["GroupLR:TimepointT2"]]
    full <- beta(d)
    shift <- vapply(unique(d$subject), function(s) full - beta(d[d$subject != s, ]), numeric(1))
    tibble(
      measure = MEASURES$label[i],
      max_abs_shift = max(abs(shift)),
      dfbetas_guide = 2 / sqrt(dplyr::n_distinct(d$subject))
    )
  }))
}
