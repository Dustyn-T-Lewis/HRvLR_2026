# Data-driven version (preferred): WGCNA module eigengenes from F06, already
# per-sample, run through the same HLM and related to the response phenotypes.
# Cleaner than curated-pathway GSVA - 11 data-driven modules, light correction.
if (!exists("expr")) source(here::here("05_test", "a_script", "setup.R"))
pacman::p_load(lme4, lmerTest)

eig <- read_csv(
  here("04_Figures", "F06", "c_data", "wgcna_eigengene.csv"),
  show_col_types = FALSE
) |>
  mutate(
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3")),
    group = factor(group_arm, levels = c("LR", "HR")),
    module = group_id
  )
mods <- unique(eig$module)

# (a) module trajectory divergence: ME ~ group * timepoint + (1|subject)
traj <- bind_rows(lapply(mods, function(m) {
  d <- eig[eig$module == m, ]
  a <- anova(lmerTest::lmer(ME ~ group * timepoint + (1 | subject), data = d))
  eff <- d |>
    group_by(group, timepoint) |>
    summarise(mu = mean(ME), .groups = "drop") |>
    pivot_wider(names_from = c(group, timepoint), values_from = mu)
  tibble(
    module = m, p_interaction = a["group:timepoint", "Pr(>F)"],
    slope = (eff$HR_T3 - eff$HR_T1) - (eff$LR_T3 - eff$LR_T1)
  )
}))
traj$fdr <- p.adjust(traj$p_interaction, "BH")
write_csv(traj, file.path(DAT, "module_trajectory.csv"))

fig_traj <- ggplot(
  traj |> mutate(
    module = reorder(module, -log10(p_interaction)), hit = fdr < 0.05
  ),
  aes(-log10(p_interaction), module)
) +
  geom_segment(aes(x = 0, xend = -log10(p_interaction), yend = module, color = module),
    linewidth = 0.6
  ) +
  geom_point(aes(fill = module, shape = hit), size = 4.5, stroke = 0.4, color = "grey30") +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "grey55") +
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
  labs(
    x = "-log10 p (group x timepoint)", y = "WGCNA module",
    title = "Module trajectory divergence (HR vs LR)",
    subtitle = sprintf("%d modules; nominal p<0.05: %d, FDR<0.05: %d", nrow(traj), sum(traj$p_interaction < 0.05), sum(traj$fdr < 0.05))
  ) +
  FIG_THEME
save_panel(fig_traj, file.path(RPT, "fig7_module_trajectory"), 150, 120)

# (b) how each module's response relates to the phenotype: eigengene change
# (T1 -> T3) per subject vs each response trait
TRAITS <- c("comp_hypertrophy", "d_mcsa", "d_fcsa_I", "d_fcsa_II", "d_1rm_legpress", "d_1rm_ext")
chg <- eig |>
  filter(timepoint %in% c("T1", "T3")) |>
  select(subject, module, timepoint, ME) |>
  pivot_wider(names_from = timepoint, values_from = ME) |>
  mutate(delta = T3 - T1)
assoc <- bind_rows(lapply(TRAITS, function(tr) {
  bind_rows(lapply(mods, function(m) {
    d <- chg[chg$module == m, ] |>
      left_join(pheno |> select(subject, all_of(tr)), by = "subject")
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
write_csv(assoc, file.path(DAT, "module_phenotype.csv"))

trait_lab <- c(
  comp_hypertrophy = "Comp. hypertrophy", d_mcsa = "d mCSA", d_fcsa_I = "d fCSA I",
  d_fcsa_II = "d fCSA II", d_1rm_legpress = "d 1RM leg", d_1rm_ext = "d 1RM ext"
)
fig_assoc <- ggplot(
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
    x = NULL, y = "WGCNA module",
    title = "How each module's response relates to the phenotype",
    subtitle = "Eigengene change (T1->T3) vs response trait; r shown if p<0.05, * = FDR<0.05"
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_panel(fig_assoc, file.path(RPT, "fig8_module_phenotype"), 165, 135)

cat(sprintf(
  "module HLM: trajectory FDR<0.05 %d; phenotype assoc p<0.05 %d\n",
  sum(traj$fdr < 0.05), sum(assoc$p < 0.05)
))
