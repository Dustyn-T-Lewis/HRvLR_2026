# Per-person views of a feature-by-biopsy matrix, and the phenotype as the trait steps read it.
# Sourced by 02_HRvLR steps 06 and 07 and 03_Continuous steps 01 and 03, after dplyr and tidyr load.

# One column per person: the value at time point `at`, or with `from`, the change from `from` to
# `at`. Every person has all three biopsies, so no one drops out.
person_values <- function(values, meta, at, from = NULL) {
  at_time <- \(tp) {
    rows <- meta$timepoint == tp
    m <- values[, meta$sample_id[rows], drop = FALSE]
    colnames(m) <- meta$subject[rows]
    m
  }
  if (is.null(from)) {
    return(at_time(at))
  }
  before <- at_time(from)
  at_time(at)[, colnames(before), drop = FALSE] - before
}

# The phenotype as six windows, each a table with one row per person: the size and strength traits
# as % change (with the composite hypertrophy score and total volume load) beside the change and T1
# windows, and as levels beside the level windows. T3 has no phenotype, so the T3 level is paired
# with the T2 traits.
trait_windows <- function(phenotype) {
  traits <- c("fcsa_mixed", "fcsa_type1", "fcsa_type2", "mcsa", "leg_press_1rm", "leg_ext_1rm")
  levels <- phenotype |>
    select(subject, matches("_(pre|post)_")) |>
    pivot_longer(-subject, names_to = c("trait", "when"), names_pattern = "(.*)_(pre|post)_.*") |>
    filter(trait %in% traits) |>
    pivot_wider(names_from = trait)
  change <- levels |>
    pivot_longer(all_of(traits), names_to = "trait") |>
    pivot_wider(names_from = when) |>
    transmute(subject, trait = paste0("pct_", trait), pct = 100 * (post - pre) / pre) |>
    pivot_wider(names_from = trait, values_from = pct) |>
    left_join(select(phenotype, subject, comp_hypertrophy, volume_load_total_kg), by = "subject")
  level_at <- \(when) select(filter(levels, .data$when == .env$when), subject, all_of(traits))
  list(
    baseline_level = level_at("pre"), baseline = change, training = change,
    post_level = level_at("post"), acute = change, post_bout_level = level_at("post")
  )
}

# The feature value each window pairs with its traits.
window_values <- function(values, meta) {
  list(
    baseline_level = person_values(values, meta, "T1"),
    baseline = person_values(values, meta, "T1"),
    training = person_values(values, meta, "T2", from = "T1"),
    post_level = person_values(values, meta, "T2"),
    acute = person_values(values, meta, "T3", from = "T2"),
    post_bout_level = person_values(values, meta, "T3")
  )
}
