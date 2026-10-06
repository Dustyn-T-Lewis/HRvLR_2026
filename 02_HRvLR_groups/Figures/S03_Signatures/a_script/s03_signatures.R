# Supplement to Figure 1. a, b: the training and acute signature ranks repeated with the smaller
# Pi < 0.05 sets. c: each arm's p < 0.05 training sets followed through T1, T2 and T3 in both
# arms, which adds the bout step that a single-window ranking cannot show. Reads saved outputs
# only.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source(here("R", "helpers.R"))
source(here("R", "figure_theme.R"))
source(here("R", "figure_panels.R"))

stage <- here("02_HRvLR_groups", "Figures", "S03_Signatures")
figure_dir <- file.path(stage, "b_reports")
out <- file.path(stage, "c_data")
walk(c(figure_dir, out), dir.create, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  contrasts = "02_HRvLR_groups/02_DEP/c_data/02_contrasts.xlsx",
  set_tests = "02_HRvLR_groups/03_Pathways/c_data/set_tests.rds"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Missing inputs: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
read_sheet <- \(sheet) readxl::read_excel(paths[["contrasts"]], sheet)
proteins <- readRDS(paths[["set_tests"]])$protein_results
fry_rows <- read_sheet("fry_concordance")

# Panels a and b.
barcodes <- map(
  set_names(c("Training", "Acute")),
  \(window) signature_barcode(proteins, fry_rows, window, set_level = "pi")
)

# Panel c: mean change from each person's T1 (02_DEP, set_trajectory). In the arm that chose
# a set, the training window is where it was selected, so it carries no test, and the bout window
# shares T2 with the selection: a protein high at T2 by chance reads as up, then down. Those cells
# are faded and the bout result marked with a dagger. Each window reports the share of proteins
# moving as they did in training and fry's FDR on that window's contrast.
trajectory <- read_sheet("set_trajectory") |>
  filter(startsWith(set, "nominal_")) |>
  mutate(
    chosen = sub(".*_", "", set_from), direction = sub(".*_", "", set), own = chosen == arm
  )
window_stats <- fry_rows |>
  filter(startsWith(set_from, "Training"), startsWith(set, "nominal_")) |>
  transmute(
    chosen = sub(".*_", "", set_from), arm = sub(".*_", "", ranked_on),
    window = sub("_.*", "", ranked_on), direction = sub(".*_", "", set),
    n_proteins, same_sign_share, fdr,
    label = sprintf("%.0f%% same way\nFDR %.2g", 100 * same_sign_share, fdr)
  )
selection_windows <- distinct(trajectory, chosen, arm, direction) |>
  filter(chosen == arm) |>
  mutate(window = "Training", label = "selection\nwindow")
window_labels <- bind_rows(window_stats, selection_windows) |>
  mutate(
    label = if_else(window == "Acute" & chosen == arm, paste0(label, " \u2020"), label),
    x = if_else(window == "Training", 1.5, 2.5),
    y = if_else(direction == "up", Inf, -Inf), vjust = if_else(direction == "up", 1.4, -0.5)
  )
cell_names <- c(
  chosen_HR = "Changed by training in HR", chosen_LR = "Changed by training in LR",
  in_HR = "Measured in HR", in_LR = "Measured in LR"
)
as_cells <- \(d) mutate(d, chosen = paste0("chosen_", chosen), arm = paste0("in_", arm))
shade_window <- \(from, to, family) {
  annotate(
    "rect",
    xmin = from, xmax = to, ymin = -Inf, ymax = Inf,
    fill = contrast_tone("fill")[[family]], alpha = 0.6
  )
}
panel_trajectory <- ggplot(
  as_cells(trajectory), aes(timepoint, mean_change, group = direction, colour = direction)
) +
  shade_window(1, 2, "Training") +
  shade_window(2, 3, "Acute") +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = mean_change - se, ymax = mean_change + se, alpha = own),
    width = 0.06, linewidth = 0.3
  ) +
  geom_line(aes(alpha = own), linewidth = 0.6) +
  geom_point(aes(alpha = own), size = 1.3) +
  geom_text(
    aes(x, y, label = label, vjust = vjust, colour = direction), as_cells(window_labels),
    inherit.aes = FALSE, size = 1.8, lineheight = 0.9
  ) +
  scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 1), guide = "none") +
  scale_colour_manual(
    values = c(up = direction_colours[["Up"]], down = direction_colours[["Down"]]),
    labels = c(down = "decreased in training", up = "increased in training"), name = NULL
  ) +
  scale_x_discrete(expand = expansion(add = 0.04)) +
  scale_y_continuous(expand = expansion(mult = 0.45)) +
  ggh4x::facet_grid2(
    chosen ~ arm,
    scales = "free_y", independent = "y", switch = "y",
    labeller = as_labeller(cell_names)
  ) +
  labs(
    x = NULL, y = NULL, title = "Where training-responsive proteins go after the bout",
    subtitle = paste(
      "p < 0.05 sets; change from own T1 (log2) \u00b1 SE; % same way, fry FDR;",
      "\u2020 T2 shared"
    )
  ) +
  figure_theme +
  theme(
    legend.position = "bottom", strip.placement = "outside",
    strip.text = element_text(size = 6), panel.spacing = unit(2, "mm")
  )

# Every plotted number traces to its source table before anything is written.
pi_numbers <- map(barcodes, "labels") |> list_rbind()
stopifnot(
  nrow(pi_numbers) == 8,
  nrow(trajectory) == 4 * 2 * 3,
  # Other arm's training window, plus both arms' bout windows, per set.
  nrow(window_stats) == 12,
  all(count(window_labels, chosen, arm, direction)$n == 2)
)

figure <- wrap_plots(
  barcodes$Training$plot, barcodes$Acute$plot, panel_trajectory,
  nrow = 1
) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, figure_dir, "s03_signatures", width = 300, height = 110)

sheets <- list(
  pi_barcode = select(pi_numbers, window, chosen, ranked_in, direction, n, same_sign_share, fdr),
  trajectory = select(trajectory, set_from, set, n_proteins, arm, timepoint, mean_change, se),
  window_stats = select(
    window_stats, chosen, arm, window, direction, n_proteins, same_sign_share, fdr
  )
)
holds <- c(
  "Panels a, b: each arm's Pi < 0.05 sets in the other arm's ranking: n, % same way, fry FDR.",
  "Panel c: mean change from each subject's T1 (and SE) per training set, arm and time point.",
  "Panel c: per cell and window, share of proteins moving as in training and fry FDR."
)
write_workbook(sheets, holds, file.path(out, "s03_signatures.xlsx"), inputs, paths)
message("wrote s03_signatures.pdf and s03_signatures.xlsx")
