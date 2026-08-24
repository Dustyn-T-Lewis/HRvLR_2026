# Panel B: the two cells that produced a survivor, against what a random split
# of the same subjects produces. The null is mostly zero with a long tail - a
# random relabelling occasionally clears BH for dozens of proteins - so the
# counts are binned rather than plotted raw, which would hide everything behind
# the bar at zero.
if (!exists("perm_null")) {
  source(here::here("04_Figures", "F03_proteome", "a_script", "setup.R"))
}
pacman::p_load(ggplot2)

HIT_BINS <- c(0, 1, 2, 5, 10, Inf)
HIT_LABELS <- c("0", "1", "2", "3-5", "6-10", ">10")

bin_hits <- function(n) {
  cut(n, breaks = c(-1, HIT_BINS), labels = HIT_LABELS, right = TRUE)
}

build_confirm <- function(perm_null, confirmation, tag = "B") {
  hit <- confirmation |> filter(.data$n_bh >= 1)
  cell_name <- function(label, contrast) {
    paste0(LABEL_NAMES[label], " · ", contrast)
  }

  nd <- perm_null |>
    inner_join(
      hit |> dplyr::select(label, contrast),
      by = c("label", "contrast")
    ) |>
    mutate(
      cell = cell_name(.data$label, .data$contrast),
      bin = bin_hits(.data$n_bh)
    ) |>
    summarise(n_perm = dplyr::n(), .by = c(cell, bin))

  obs <- hit |>
    mutate(
      cell = cell_name(.data$label, .data$contrast),
      bin = bin_hits(.data$n_bh)
    )

  p <- ggplot(nd, aes(.data$bin, .data$n_perm)) +
    geom_col(fill = "grey80", width = 0.8) +
    geom_col(
      data = inner_join(
        nd, obs |> dplyr::select(cell, bin),
        by = c("cell", "bin")
      ),
      fill = "#B2182B", width = 0.8
    ) +
    geom_text(
      data = obs,
      aes(x = Inf, y = Inf, label = sprintf("p = %.3f", .data$p_count)),
      colour = "#B2182B", size = 2.6, hjust = 1.15, vjust = 1.6
    ) +
    facet_wrap(~cell) +
    labs(
      title = "The survivors against a shuffled-label null", tag = tag,
      subtitle = paste(
        "999 random splits of the same subjects; red = the bin the observed",
        "count falls in"
      ),
      x = "Proteins at BH < 0.05", y = "Permutations"
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      panel.grid.major.x = element_blank()
    )

  list(plot = p, audit = confirmation)
}

b_confirm <- build_confirm(perm_null, confirmation)
save_panel(
  b_confirm$plot, file.path(F03_RPT, "panels", "panel_b_confirm"),
  150, 90
)
F03_PANELS[["confirm"]] <- b_confirm$plot
F03_AUDIT[["permutation_confirmation"]] <- b_confirm$audit
