# Per-person views of a feature-by-biopsy matrix, and the phenotype as the trait steps read it.
# Expects dplyr, tidyr and purrr attached.

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

# The three feature types the trait and classification steps test, as feature-by-biopsy matrices
# in biopsy order, with each feature's readable label (gene, set name, module colour) and each set's
# collection.
load_features <- function(corrected, scores, catalog, modules, meta) {
  features <- list(
    proteins = corrected$E, sets = scores[catalog$set_id, ], modules = modules$eigengenes
  )
  walk(features, \(x) stopifnot(identical(colnames(x), meta$sample_id)))
  list(
    features = features,
    label_of = c(
      set_names(sub(";.*", "", corrected$genes$Genes), rownames(corrected)),
      set_names(catalog$pathway, catalog$set_id), set_names(rownames(modules$eigengenes))
    ),
    collection_of = set_names(catalog$collection, catalog$set_id)
  )
}
