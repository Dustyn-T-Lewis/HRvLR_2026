# F06 Panel C: how each module's response relates to the phenotype. Module
# eigengene change (T1 -> T3) per subject vs each response trait; the data-driven
# read of "which co-expression programs track the size of the response".
if (!exists("wgcna_eig")) {
  source(here::here("04_Figures", "F06", "a_script", "HRvLR_F06_setup.R"))
}

TRAITS <- c("comp_hypertrophy", "d_mcsa", "d_fcsa_I", "d_fcsa_II", "d_1rm_legpress", "d_1rm_ext")
mods <- setdiff(unique(wgcna_eig$group_id), "grey")

chg <- wgcna_eig |>
  dplyr::filter(group_id %in% mods, timepoint %in% c("T1", "T3")) |>
  dplyr::select(subject, module = group_id, timepoint, ME) |>
  tidyr::pivot_wider(names_from = timepoint, values_from = ME) |>
  dplyr::mutate(delta = T3 - T1)

assoc <- bind_rows(lapply(TRAITS, function(tr) {
  bind_rows(lapply(mods, function(m) {
    d <- chg |>
      dplyr::filter(module == m) |>
      dplyr::left_join(pheno_tbl |> dplyr::select(subject, dplyr::all_of(tr)), by = "subject")
    d <- d[!is.na(d$delta) & !is.na(d[[tr]]), ]
    if (nrow(d) < 6) {
      return(NULL)
    }
    ct <- suppressWarnings(cor.test(d$delta, d[[tr]], method = "spearman"))
    tibble(module = m, trait = tr, r = unname(ct$estimate), p = ct$p.value, n = nrow(d))
  }))
})) |>
  group_by(trait) |>
  mutate(fdr = p.adjust(p, "BH")) |>
  ungroup()
write_csv(assoc, file.path(DAT_DIR, "module_phenotype.csv"))

panel_c <- ggplot(
  assoc |> mutate(trait = factor(trait_lab[trait], levels = trait_lab)),
  aes(trait, module, fill = r)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(p < 0.05, sprintf("%.2f%s", r, ifelse(fdr < 0.05, "*", "")), "")),
    size = 2.4, color = "grey10"
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-1, 1), name = "Spearman r"
  ) +
  labs(
    x = NULL, y = "WGCNA module", tag = "C",
    title = "Module response vs phenotype",
    subtitle = "Eigengene change (T1->T3) vs trait; r if p<0.05, * = FDR<0.05"
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)
save_panel(panel_c, file.path(RPT_DIR, "supp", "panel_c_module_phenotype"), 170, 140)
