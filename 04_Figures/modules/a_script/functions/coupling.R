pacman::p_load(dplyr, tidyr, tibble, purrr)

.phase_eigengene <- function(eigengenes, sample_meta, phase) {
  tps <- switch(phase,
    baseline = "T1",
    training = c("T1", "T2"),
    acute = c("T2", "T3")
  )
  me <- as.data.frame(eigengenes)
  me$sample_id <- rownames(me)
  long <- me |>
    pivot_longer(-sample_id, names_to = "module", values_to = "me") |>
    left_join(sample_meta, by = "sample_id") |>
    filter(timepoint %in% tps)
  if (phase == "baseline") {
    long |> select(subject, module, value = me)
  } else {
    long |>
      pivot_wider(id_cols = c(subject, module), names_from = timepoint, values_from = me) |>
      mutate(value = .data[[tps[2]]] - .data[[tps[1]]]) |>
      select(subject, module, value)
  }
}

module_coupling <- function(eigengenes, sample_meta, pheno, phase) {
  traits <- c("comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa", "d_1rm_legpress", "d_1rm_ext")
  traits <- intersect(traits, names(pheno))
  me_phase <- .phase_eigengene(eigengenes, sample_meta, phase)

  grid <- tidyr::expand_grid(module = unique(me_phase$module), trait = traits)
  res <- purrr::pmap_dfr(grid, function(module, trait) {
    d <- me_phase |>
      filter(module == !!module) |>
      left_join(pheno |> select(subject, all_of(trait)), by = "subject")
    ct <- suppressWarnings(cor.test(d$value, d[[trait]], method = "pearson"))
    tibble(module = module, trait = trait, r = unname(ct$estimate), p = ct$p.value)
  })
  res |> mutate(p_bh = p.adjust(p, method = "BH"), phase = phase)
}
