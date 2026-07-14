# Panel D: does a module's baseline eigengene predict how much a subject adapts?
#
# It does not, and the panel says so on its face. Two numbers are reported side by side for
# every module x trait cell:
#
#   r      the in-sample correlation, which is what a naive reading would quote
#   q2     leave-one-subject-out cross-validated R2; q2 <= 0 means the model predicts the
#          held-out subject WORSE than the trait's own mean does
#
# The gap between them is the point. An earlier version of this analysis reported a q2 of
# 0.713 - an artifact of modules that encoded subject identity (see wgcna.R). With the
# modules defined on within-subject co-regulation, the best cell in the whole 12 x 6 grid is
# a q2 of ~0.08, and that is the ARGMAX of 72 cells, so it is nearly free.
#
# This panel owns its supplement: the mixed model that keeps the full repeated-measures
# structure this panel collapses to T1.
if (!exists("mods")) source(here::here("04_Figures", "F04_modules", "a_script", "setup.R"))
pacman::p_load(ggplot2, dplyr, tidyr, purrr, forcats, lmerTest, lme4)

loo_q2 <- function(x, y) {
  d <- data.frame(x = x, y = y)
  d <- d[stats::complete.cases(d), ]
  if (nrow(d) < 6) {
    return(NA_real_)
  }
  pred <- vapply(
    seq_len(nrow(d)),
    function(i) stats::predict(stats::lm(y ~ x, d[-i, ]), d[i, ]),
    numeric(1)
  )
  1 - sum((d$y - pred)^2) / sum((d$y - mean(d$y))^2)
}

t1 <- me_long |>
  filter(timepoint == "T1") |>
  select(subject, module, ME)

grid <- tidyr::expand_grid(module = unique(t1$module), trait = ADAPTATION_TRAITS)

pred <- purrr::pmap_dfr(grid, function(module, trait) {
  d <- t1 |>
    filter(module == !!module) |>
    left_join(pheno |> select(subject, value = all_of(trait)), by = "subject") |>
    filter(!is.na(ME), !is.na(value))
  ct <- suppressWarnings(stats::cor.test(d$ME, d$value))
  tibble(
    module = module, trait = trait, n = nrow(d),
    r = unname(ct$estimate), p = ct$p.value, q2 = loo_q2(d$ME, d$value)
  )
}) |>
  mutate(fdr = p.adjust(p, "BH"))

pred_plot <- pred |>
  mutate(
    module = sub("^ME", "", module),
    module = fct_reorder(module, q2, .fun = max),
    trait = factor(trait, levels = ADAPTATION_TRAITS)
  )

pD <- ggplot(pred_plot, aes(trait, module, fill = q2)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", q2)), size = 2.1, color = "grey15") +
  scale_fill_gradient2(
    low = "#B2182B", mid = "grey95", high = "#2166AC", midpoint = 0,
    name = expression(q^2), limits = c(-1, 1), oob = scales::squish
  ) +
  scale_x_discrete(labels = TRAIT_LABELS) +
  labs(
    title = "Can a baseline module predict the size of the response?",
    subtitle = sprintf(
      "Leave-one-subject-out q2. Blue = predicts; red = worse than the trait's own mean. Best of %d cells: q2 = %.3f. Nothing survives FDR (min q = %.2f).",
      nrow(pred), max(pred$q2, na.rm = TRUE), min(pred$fdr, na.rm = TRUE)
    ),
    x = NULL, y = NULL, tag = "D"
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40")
  )

save_panel(pD, file.path(RPT_DIR, "panels", "panel_d_phenotype"), 190, 120)
F04_PANELS[["phenotype"]] <- pD
F04_AUDIT[["module_prediction"]] <- pred

# Supplement: the mixed model. Panel D collapses the design to T1; this keeps every
# timepoint, fitting ME ~ trait * timepoint + (1 | subject). The trait:timepoint term asks
# whether a module's whole time-course bends with the size of the response - the question
# Panel D can only ask at one point in time. Singular fits are flagged, not dropped: a
# time-invariant trait can drive the random-intercept variance to zero, which inflates the
# trait main effect, so the interaction is the term to trust.
fit_cell <- function(module, trait) {
  d <- me_long |>
    filter(module == !!module) |>
    left_join(pheno |> select(subject, trait_value = all_of(trait)), by = "subject") |>
    filter(!is.na(ME), !is.na(trait_value))
  fit <- tryCatch(
    suppressWarnings(
      lmerTest::lmer(ME ~ trait_value * timepoint + (1 | subject), data = d)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(
      module = module, trait = trait, n = nrow(d),
      p_trait = NA_real_, p_trait_time = NA_real_, singular = NA
    ))
  }
  a <- stats::anova(fit)
  tibble(
    module = module, trait = trait, n = nrow(d),
    p_trait = a["trait_value", "Pr(>F)"],
    p_trait_time = a["trait_value:timepoint", "Pr(>F)"],
    singular = lme4::isSingular(fit)
  )
}

lmm <- purrr::pmap_dfr(grid, fit_cell) |>
  group_by(trait) |>
  mutate(
    fdr_trait = p.adjust(p_trait, "BH"),
    fdr_trait_time = p.adjust(p_trait_time, "BH")
  ) |>
  ungroup()

lmm_plot <- lmm |>
  mutate(
    module = sub("^ME", "", module),
    trait = factor(trait, levels = ADAPTATION_TRAITS)
  )

p_lmm <- ggplot(lmm_plot, aes(trait, module, fill = -log10(fdr_trait_time))) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "grey95", high = "#6A51A3", name = expression(-log[10] * " BH q")) +
  scale_x_discrete(labels = TRAIT_LABELS) +
  labs(
    title = "Does the module's time-course bend with the response? (trait x timepoint)",
    subtitle = sprintf(
      "ME ~ trait * timepoint + (1 | subject), BH within trait. Minimum q = %.2f. Nothing survives.",
      min(lmm$fdr_trait_time, na.rm = TRUE)
    ),
    x = NULL, y = NULL
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40")
  )

save_panel(p_lmm, file.path(RPT_DIR, "supp", "panel_d_supp_phenotype_lmm"), 190, 120)
F04_AUDIT[["module_phenotype_lmm"]] <- lmm
