# F02 Panel G: supervised separation at baseline (sPLS-DA)
# The supervised counterpart to Panel A. sPLS-DA is handed the arm labels and asked
# for the axis that best splits them, so at 1898 proteins and 15 subjects it separates
# whatever it is given. The right-hand facet fits the same pipeline to labels shuffled
# across subjects and separates just as cleanly; the held-out numbers come from the
# leave-one-subject-out screen, not from this fit. Title on the composite.

pacman::p_load(here, dplyr, tibble, ggplot2, mixOmics, readr, withr)

if (!exists("meta")) source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))

PG_W <- 120
PG_H <- 110
PG_KEEPX <- c(50, 50)
PG_SEED <- 42

t1 <- meta |> filter(Timepoint == "T1", Col_ID %in% colnames(imp_mat))
x_t1 <- t(imp_mat[, t1$Col_ID])
x_t1 <- x_t1[, apply(x_t1, 2, var) > 0]

# Shuffle the arm label once per subject. Permuting sample rows would break the
# within-subject structure and give a null that is too easy to beat.
shuffled <- with_seed(PG_SEED, {
  subj <- unique(t1$subject)
  stats::setNames(sample(t1$Group[match(subj, t1$subject)]), subj)
})

splsda_scores <- function(y, facet) {
  fit <- mixOmics::splsda(x_t1, y, ncomp = 2, keepX = PG_KEEPX)
  as.data.frame(fit$variates$X) |>
    rlang::set_names(c("comp1", "comp2")) |>
    mutate(fitted_label = y, facet = facet)
}

scores <- bind_rows(
  splsda_scores(t1$Group, "Arm labels"),
  splsda_scores(unname(shuffled[t1$subject]), "Labels shuffled")
) |>
  mutate(facet = factor(facet, levels = c("Arm labels", "Labels shuffled")))

# Held-out performance for this exact cell, carried from the screen store so the
# panel cannot drift from the number that qualifies it.
cv <- read_csv(
  here("04_Figures", "supp_screens", "c_data", "class_cells_summary.csv"),
  show_col_types = FALSE
) |>
  filter(
    .data$B == 200, .data$level == "proteins",
    .data$config == "T1", .data$model == "spls"
  )

# Only the left facet has a held-out number to report. Repeating it under the
# shuffle would imply the shuffle was cross-validated too.
notes <- tibble(
  facet = factor(levels(scores$facet), levels = levels(scores$facet)),
  note = c(
    sprintf(
      "Leave-one-subject-out\nAUC = %.2f (p = %.2f)",
      cv$estimate[1], cv$perm_p[1]
    ),
    "Same pipeline,\narm label shuffled"
  )
)

pG <- ggplot(scores, aes(comp1, comp2, colour = fitted_label, fill = fitted_label)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.8, alpha = 0.15, linewidth = 0.3) +
  geom_point(size = 2.1) +
  facet_wrap(~facet) +
  scale_colour_manual(values = GROUP_COLORS) +
  scale_fill_manual(values = GROUP_COLORS) +
  geom_label(
    data = notes, aes(x = -Inf, y = Inf, label = note),
    hjust = -0.04, vjust = 1.08, inherit.aes = FALSE,
    size = FIG_GEOM_TEXT - 0.7, fontface = "bold", colour = "grey15",
    fill = scales::alpha("white", 0.85),
    label.size = 0, label.padding = unit(2, "pt"), lineheight = 0.95
  ) +
  labs(x = "sPLS-DA component 1", y = "sPLS-DA component 2", colour = NULL, fill = NULL) +
  FIG_THEME +
  theme(legend.position = "bottom", plot.margin = margin(t = 6, r = 5, b = 4, l = 4))

save_png(pG, file.path(RPT_DIR, "panels", "panel_g_splsda"), PG_W, PG_H)
F02_AUDIT[["panel_G_splsda_scores"]] <- scores |>
  mutate(subject = rep(t1$subject, 2), true_group = rep(t1$Group, 2))
cat("F02 Panel G done.\n")
