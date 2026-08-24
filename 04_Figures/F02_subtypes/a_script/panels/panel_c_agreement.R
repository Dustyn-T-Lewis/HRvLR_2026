# Panel C: if a two-group split is forced anyway, does it land on any label the
# study already has? Labelled descriptive because the gate in panel A stayed
# shut; it would be the headline if it had opened.
if (!exists("forced_two_group")) {
  source(here::here("04_Figures", "F02_subtypes", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats)

LABEL_NAMES <- c(
  given = "Given HR/LR", fcsa_I = "fCSA type I", fcsa_II = "fCSA type II",
  mcsa = "Whole-muscle CSA", `1rm_legpress` = "1RM leg press",
  `1rm_ext` = "1RM leg extension"
)

build_agreement <- function(forced_two_group, tag = "C") {
  d <- forced_two_group |>
    mutate(
      space_lab = SPACE_LABELS[.data$space],
      label_lab = fct_rev(factor(LABEL_NAMES[.data$label],
        levels = unname(LABEL_NAMES)
      ))
    )

  p <- ggplot(d, aes(.data$ari, .data$label_lab, fill = .data$space_lab)) +
    geom_vline(
      xintercept = 0, linetype = "dashed", colour = "grey50",
      linewidth = 0.4
    ) +
    geom_col(position = position_dodge(0.75), width = 0.65) +
    scale_fill_manual(values = unname(GROUP_COLORS), name = NULL) +
    scale_x_continuous(limits = c(-0.25, 1)) +
    labs(
      title = "A forced two-group split matches nothing", tag = tag,
      subtitle = paste(
        "Adjusted Rand index vs each label;",
        "0 = chance agreement, 1 = identical"
      ),
      x = "Adjusted Rand index", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "top",
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = d |> dplyr::select(space, label, n, ari, gate_open))
}

c_agreement <- build_agreement(forced_two_group)
save_panel(
  c_agreement$plot, file.path(F02_RPT, "panels", "panel_c_agreement"),
  160, 95
)
F02_PANELS[["agreement"]] <- c_agreement$plot
F02_AUDIT[["forced_two_group_agreement"]] <- c_agreement$audit
