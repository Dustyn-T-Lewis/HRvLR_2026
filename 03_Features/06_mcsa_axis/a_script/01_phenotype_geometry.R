#!/usr/bin/env Rscript
# Is the six-trait phenotype one axis or two, and how much of the answer is one
# subject?
#
# Three readouts on the 16 subjects: the pairwise Spearman structure, the same
# structure with LR_S14 removed, and the composite's own decomposition onto the
# three CSA measures. The decomposition matters because it settles what the
# HR/LR split was made of: d_mcsa is a weighted ingredient of the composite,
# not a phenotype outside it, so a d_mcsa result cannot be read as answering a
# different question than F04 asked.
#
# Read 01_phenotype_geometry.csv for the correlations both ways and
# 01_composite_weights.csv for the formula weights beside the galamm loadings
# the measurement model fitted. A discordant subject looks like a correlation
# that moves a long way when one row leaves.

pacman::p_load(here, dplyr, tidyr, readr, tibble, purrr)
source(here("03_Features", "06_mcsa_axis", "a_script", "mcsa_helpers.R"))

OUT <- here("03_Features", "06_mcsa_axis", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

pheno <- phenotype_table()

pairwise_rho <- function(d, subject_set) {
  m <- stats::cor(
    d[, TRAITS], method = "spearman", use = "pairwise.complete.obs"
  )
  as.data.frame(as.table(m)) |>
    rlang::set_names(c("trait_x", "trait_y", "rho")) |>
    filter(as.integer(.data$trait_x) < as.integer(.data$trait_y)) |>
    mutate(subject_set = subject_set)
}

geometry <- bind_rows(
  pairwise_rho(pheno, "all 16"),
  pairwise_rho(
    filter(pheno, .data$subject != DISCORDANT_SUBJECT), "LR_S14 dropped"
  )
) |>
  pivot_wider(names_from = "subject_set", values_from = "rho") |>
  mutate(shift = .data$`LR_S14 dropped` - .data$`all 16`) |>
  arrange(desc(abs(.data$shift)))

write_csv(geometry, file.path(OUT, "01_phenotype_geometry.csv"))

loadings <- read_csv(
  here("03_Features", "04_galamm_pilot", "c_data", "02_q2_measurement.csv"),
  show_col_types = FALSE
)

weights <- composite_weights(pheno) |>
  left_join(
    select(loadings, "item", "loading", "var_explained"),
    by = "item"
  ) |>
  mutate(weight_share = .data$weight / sum(.data$weight))

write_csv(weights, file.path(OUT, "01_composite_weights.csv"))

separation <- tibble(
  trait = TRAITS,
  auc = map_dbl(TRAITS, \(t) arm_auc(pheno[[t]], pheno$group_arm)),
  wilcox_p = map_dbl(TRAITS, \(t) {
    stats::wilcox.test(pheno[[t]] ~ pheno$group_arm)$p.value
  })
)
write_csv(separation, file.path(OUT, "01_arm_separation.csv"))

trimmed <- filter(pheno, .data$subject != DISCORDANT_SUBJECT)
cat(sprintf(
  paste0(
    "composite R2 %.3f with mCSA, %.3f without; mCSA weight %.2f of %.2f\n",
    "d_mcsa vs fibre axis: rho %.2f all 16, %.2f without %s\n",
    "d_mcsa alone separates the arms at AUC %.2f (wilcoxon p = %.3f)\n"
  ),
  weights$r2_full[1], weights$r2_without_mcsa[1],
  weights$weight[weights$item == "d_mcsa"], sum(weights$weight),
  stats::cor(pheno$d_mcsa, fibre_axis(pheno), method = "spearman"),
  stats::cor(trimmed$d_mcsa, fibre_axis(trimmed), method = "spearman"),
  DISCORDANT_SUBJECT,
  separation$auc[separation$trait == "d_mcsa"],
  separation$wilcox_p[separation$trait == "d_mcsa"]
))
