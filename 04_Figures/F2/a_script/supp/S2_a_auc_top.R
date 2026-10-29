library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))
source(here("04_Figures", "style", "module_labels.R"))

inputs <- c(
  classification = "02_HRvLR/06_Classification/c_data/classification.rds",
  characterisation = characterisation_input
)
classification <- readRDS(here(inputs[["classification"]]))
modules <- module_identity(read_module_labels(), width = 42)

# The five tasks sit under the contrast strips: the levels under Between groups, the changes under
# their window.
tasks <- c(
  Baseline_HRvLR = "T1", Trained_HRvLR = "T2", Acute_HRvLR = "T3",
  Training_change_HRvLR = "T2 - T1", Acute_change_HRvLR = "T3 - T2"
)
task_family <- c(
  Baseline_HRvLR = "Between groups", Trained_HRvLR = "Between groups",
  Acute_HRvLR = "Between groups", Training_change_HRvLR = "Training",
  Acute_change_HRvLR = "Acute"
)

# The step's top features per task and feature type (smallest Wilcoxon p, then AUC furthest from
# 0.5).
counts <- classification$counts
n_top <- max(classification$top$rank)
panel_data <- classification$top |>
  left_join(
    select(counts, feature_type, task, n_tested, n_nominal),
    by = c("feature_type", "task")
  ) |>
  select(
    feature_type, task, rank, feature, collection, label, auc, auc_lower, auc_upper, p, fdr,
    n_nominal, n_tested
  )

check <- inner_join(panel_data, classification$auc, by = c("feature_type", "task", "feature"))
stopifnot(
  nrow(check) == nrow(panel_data), all(check$auc.x == check$auc.y),
  all(check$p.x == check$p.y), setequal(names(tasks), panel_data$task),
  setequal(names(feature_types), panel_data$feature_type),
  all(count(panel_data, feature_type, task)$n == n_top)
)

d <- panel_data |>
  mutate(
    task_row = factor(tasks[task], tasks),
    type = factor(feature_types[feature_type], feature_types),
    name = case_when(
      feature_type == "sets" ~ sprintf(
        "%s (%s)", pathway_label(label, 42), collection_names[collection]
      ),
      feature_type == "modules" ~ modules[feature],
      .default = label
    ),
    row = paste(task, feature_type, rank),
    evidence = evidence_level(p, fdr),
    higher = factor(if_else(auc > 0.5, "HR", "LR"), c("HR", "LR"))
  )
tallies <- distinct(d, task_row, type, n_nominal, n_tested) |>
  mutate(
    chance = n_tested * alpha_nominal,
    text = sprintf(
      "%d at p < %s, %s by chance", n_nominal, alpha_nominal,
      if_else(chance < 10, round(chance, 1), round(chance))
    )
  )
# Colour already gives the direction, so the evidence code here is a shape alone: filled at FDR,
# open at p only, a cross otherwise.
auc_shapes <- rlang::set_names(c(16, 1, 4), evidence_levels)

families <- unname(task_family[names(tasks)])
panel <- ggplot(d, aes(auc, reorder(row, -rank))) +
  # The chance line stops above the tally line below the bottom row.
  annotate(
    "segment",
    x = 0.5, xend = 0.5, y = 0.7, yend = Inf,
    linetype = "dashed", linewidth = 0.3, colour = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = auc_lower, xmax = auc_upper),
    orientation = "y", width = 0.3, linewidth = 0.35,
    na.rm = TRUE
  ) +
  geom_point(aes(colour = higher, shape = evidence), size = 1.5, stroke = 0.5) +
  geom_text(
    aes(x = 0, y = 0.35, label = text), tallies,
    inherit.aes = FALSE, hjust = 0, size = 1.9, colour = "grey30"
  ) +
  ggh4x::facet_grid2(
    task_row ~ type,
    scales = "free_y", independent = "y", switch = "y", strip = contrast_strips(y = families)
  ) +
  scale_y_discrete(
    labels = tibble::deframe(distinct(d, row, name)), expand = expansion(add = c(1.1, 0.6))
  ) +
  scale_colour_manual(
    values = group_colours, labels = c("higher in HR", "higher in LR"),
    name = NULL
  ) +
  scale_shape_manual(values = auc_shapes, name = NULL, drop = TRUE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Can one feature tell HR from LR?",
    subtitle = sprintf(
      paste(
        "Top %d features per task and feature type by Wilcoxon p; AUC with DeLong 95%% CI",
        "(no interval at AUC 0 or 1); dashed: chance"
      ),
      n_top
    ),
    x = "AUC (above 0.5: higher in HR)", y = NULL
  ) +
  figure_theme +
  theme(
    strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 5.5),
    strip.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 5), panel.grid.major.y = element_blank(),
    legend.position = "bottom", panel.spacing.x = unit(3, "mm"),
    panel.spacing.y = unit(1.5, "mm")
  )

save_panel(panel, "S2", "a_auc_top", 250, 105)
