# F06 supplement: module trajectory moderation by phenotype. The LMM companion
# to Panel C. Where Panel C correlates the T1->T3 eigengene delta with each trait,
# this fits ME ~ trait * timepoint + (1 | subject) and reads the trait:timepoint
# interaction, the test of whether a module's whole time-course bends with the
# size of the response. Shown as the null it is: nothing survives FDR.
if (!exists("pheno_lmm")) {
  source(here::here("04_Figures", "F06", "a_script", "HRvLR_F06_setup.R"))
}

lmm_grid <- pheno_lmm |>
  dplyr::filter(trait %in% names(trait_lab)) |>
  dplyr::mutate(
    neg_log10_p = -log10(p_trait_time),
    trait = factor(trait_lab[trait], levels = trait_lab)
  )

n_nominal <- sum(lmm_grid$p_trait_time < 0.05, na.rm = TRUE)
n_fdr <- sum(lmm_grid$fdr_trait_time < 0.05, na.rm = TRUE)

panel_s <- ggplot(lmm_grid, aes(trait, module, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(
      p_trait_time < 0.05,
      sprintf("%.3f%s", p_trait_time, ifelse(fdr_trait_time < 0.05, "*", "")),
      ""
    )),
    size = 2.3, color = "grey10"
  ) +
  scale_fill_gradient(
    low = "white", high = "#4A1486", limits = c(0, NA),
    name = expression(-log[10] ~ italic(p))
  ) +
  labs(
    x = NULL, y = "WGCNA module", tag = "S",
    title = "Module trajectory moderation by phenotype (LMM)",
    subtitle = sprintf(
      "trait x timepoint interaction; %d/%d nominal p<0.05, %d at FDR<0.05",
      n_nominal, nrow(lmm_grid), n_fdr
    )
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)
save_panel(panel_s, file.path(RPT_DIR, "supp", "panel_s_phenotype_lmm"), 170, 140)
