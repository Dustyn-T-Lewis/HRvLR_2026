# S7: module detail from 02_HRvLR/08_Module_Detail, five PDFs beside the other supplements: hubs and
# rewiring (one page per module), coupling within each group, each person's response, and the PLS
# models. It draws the step's tables whole, so it keeps no workbook of its own.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "module_labels.R"))

detail <- readRDS(here("02_HRvLR/08_Module_Detail/c_data/module_detail.rds"))
labels <- read_module_labels()
identity <- module_identity(labels, width = 45)
# One PDF per page set, every page the same size.
save_pages <- function(plots, name, width, height) {
  dir.create(figure_dir("S7", "b_reports"), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(
    file.path(figure_dir("S7", "b_reports"), paste0("S7_", name, ".pdf")),
    width = width / 25.4, height = height / 25.4, onefile = TRUE
  )
  invisible(lapply(plots, print))
  invisible(grDevices::dev.off())
}
stopifnot(setequal(unique(detail$hubs$module), labels$module))

# Hubs: each module's 15 most central proteins in the training and acute contrasts.
hub_columns <- c(
  Training_HR = "HR", Training_LR = "LR", Training_Interaction = "HR - LR",
  Acute_HR = "HR", Acute_LR = "LR", Acute_Interaction = "HR - LR"
)
hub_page <- function(m) {
  d <- filter(detail$hubs, module == m) |>
    mutate(
      gene = reorder(gene, -hub_rank),
      window = factor(sub("_.*", "", contrast), c("Training", "Acute")),
      column = factor(hub_columns[contrast], unique(hub_columns)),
      call = factor(if_else(p < alpha_nominal, "p < 0.05", "p ≥ 0.05"), c("p < 0.05", "p ≥ 0.05"))
    )
  ggplot(d, aes(column, gene)) +
    geom_point(aes(size = pmin(-log10(p), 4), fill = logFC, shape = call),
      colour = "grey30",
      stroke = 0.3
    ) +
    geom_point(
      data = filter(d, fdr < alpha_fdr), aes(size = pmin(-log10(p), 4)),
      shape = 1, colour = "black", stroke = 0.8
    ) +
    scale_fill_distiller(
      palette = "PiYG", direction = 1, limits = c(-1, 1), oob = scales::squish,
      name = "logFC\n(capped at ±1)"
    ) +
    scale_shape_manual(values = c(`p < 0.05` = 21, `p ≥ 0.05` = 1), name = NULL, drop = FALSE) +
    scale_size_area(max_size = 3, limits = c(0, 4), breaks = c(1, 2, 4), name = "-log10 p") +
    facet_grid(. ~ window) +
    labs(
      title = sprintf("Hub proteins: %s", identity[[m]]),
      subtitle = "15 highest-kME members, most central at the top; black ring: FDR < 0.05",
      x = NULL, y = NULL
    ) +
    figure_theme +
    theme(strip.text = element_text(face = "bold"))
}

# Rewiring: each member's kME within LR against within HR.
rewiring_page <- function(m) {
  d <- filter(detail$rewiring, module == m)
  named <- slice_max(d, abs(shift), n = 8)
  ggplot(d, aes(kme_LR, kme_HR)) +
    geom_abline(linetype = "22", colour = "grey50") +
    geom_point(size = 0.8, colour = "grey60") +
    geom_point(data = named, size = 1.3, colour = group_colours[["HR"]]) +
    ggrepel::geom_text_repel(
      aes(label = gene), named,
      size = 2, max.overlaps = Inf, seed = 1, segment.size = 0.2
    ) +
    coord_equal(xlim = c(-0.6, 1), ylim = c(-0.6, 1)) +
    labs(
      title = sprintf("Rewiring: %s", identity[[m]]),
      subtitle = paste(
        "Each member's kME within LR's 24 biopsies against within HR's;\non the diagonal it moves",
        "with the module equally in both; named: the 8 largest shifts"
      ),
      x = "kME in LR", y = "kME in HR"
    ) +
    figure_theme
}

save_pages(lapply(labels$module, hub_page), "hubs", 130, 100)
save_pages(lapply(labels$module, rewiring_page), "rewiring", 120, 125)

# Coupling: one lower triangle per group and window, then HR - LR boxed at Fisher p < 0.05.
module_names <- unname(identity[labels$module])
coupling <- detail$coupling |>
  pivot_longer(c(rho_HR, rho_LR, difference), names_to = "panel", values_to = "value") |>
  mutate(
    panel = factor(
      panel, c("rho_HR", "rho_LR", "difference"), c("within HR", "within LR", "HR - LR")
    ),
    window = factor(window, c("training", "acute"), c("Training change", "Acute change")),
    first = match(module_a, labels$module), second = match(module_b, labels$module),
    x = factor(module_names[pmin(first, second)], module_names),
    y = factor(module_names[pmax(first, second)], rev(module_names)),
    boxed = panel == "HR - LR" & p_difference < alpha_nominal
  )
coupling_page <- ggplot(coupling, aes(x, y, fill = value)) +
  geom_tile(colour = "white") +
  geom_tile(data = filter(coupling, boxed), fill = NA, colour = "black", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", value)), size = 1.7) +
  scale_fill_gradient2(
    low = "#4393C3", high = "#D6604D", limits = c(-1.5, 1.5), oob = scales::squish,
    name = "rho, or\nHR - LR"
  ) +
  facet_grid(window ~ panel) +
  labs(
    title = "Which modules change together, within each group",
    subtitle = sprintf(
      paste(
        "Spearman rho between two modules' eigengene changes across 8 people;",
        "box: HR - LR at Fisher p < %s"
      ),
      alpha_nominal
    ),
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank(), strip.text = element_text(face = "bold")
  )
save_pages(list(coupling_page), "coupling", 250, 200)

# Per person: every person's module changes, the people ordered by growth.
person_page <- detail$per_person |>
  mutate(
    person = reorder(sprintf("%s (%s)", subject, group), comp_hypertrophy),
    module = factor(module, labels$module, module_names),
    window = factor(window, c("training", "acute"), c("Training (T2 - T1)", "Acute (T3 - T2)"))
  ) |>
  ggplot(aes(module, person, fill = change)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#762A83", high = "#1B7837", limits = c(-0.3, 0.3), oob = scales::squish,
    name = "eigengene\nchange"
  ) +
  facet_grid(. ~ window) +
  labs(
    title = "Each person's module response",
    subtitle = "People ordered by the growth score (top grew most); group in brackets",
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 5), panel.grid = element_blank()
  )
save_pages(list(person_page), "per_person", 250, 130)

# PLS: each model's held-out score against its permutation null, then its loadings.
pls <- detail$pls |>
  mutate(task = sprintf("%s · %s · %s", model, window, outcome), task = factor(task, rev(task)))
score_plot <- ggplot(pls, aes(observed, task)) +
  geom_segment(
    aes(x = null_median, xend = null_95, yend = task),
    linewidth = 2.5, colour = "grey85"
  ) +
  geom_point(aes(colour = p_perm < 0.1), size = 2) +
  geom_text(aes(label = sprintf("p %.2f", p_perm)), nudge_y = 0.35, size = 1.8) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey55"), guide = "none") +
  ggh4x::facet_grid2(statistic ~ ., scales = "free", space = "free_y", independent = "x") +
  labs(
    title = "Do the modules together separate the groups or predict growth?",
    subtitle = paste(
      "Point: held-out score (4-fold CV repeated 5 times, balanced folds);",
      "grey bar: permutation null, median to 95th percentile"
    ),
    x = "score", y = NULL
  ) +
  figure_theme
loading_plot <- detail$pls_loadings |>
  mutate(
    task = factor(sprintf("%s · %s · %s", model, window, outcome), levels(pls$task)),
    module = factor(module, labels$module, module_names)
  ) |>
  left_join(distinct(pls, task, statistic), by = "task") |>
  ggplot(aes(module, task, fill = loading)) +
  geom_tile(colour = "white") +
  ggh4x::facet_grid2(statistic ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#762A83", high = "#1B7837", limits = c(-1, 1), name = "loading") +
  labs(x = NULL, y = NULL, subtitle = "Component-1 loadings on all 16 people") +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 5), axis.text.y = element_blank(),
    panel.grid = element_blank()
  )
save_pages(list(score_plot | loading_plot), "pls", 250, 110)
