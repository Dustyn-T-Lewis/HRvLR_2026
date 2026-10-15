# Each person's training response, trait by trait: % change in fibre area, muscle area and
# strength, the composite hypertrophy score and total volume load, as the training, T1 and acute
# windows of 03_Continuous pair them with the proteome.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("02_HRvLR", "helpers.R"))

inputs <- c(phenotype = "00_Setup/Input/phenotype.csv", helpers = "02_HRvLR/helpers.R")
phenotype <- readr::read_csv(here(inputs[["phenotype"]]), show_col_types = FALSE)
changes <- trait_windows(phenotype)$training
holds <- paste(
  "Panel a: each person's T1 to T2 change per trait (% change; composite score; volume load)."
)

traits <- setdiff(names(changes), "subject")
# The two traits on their own scales get their own facet.
own_scale <- c("comp_hypertrophy", "volume_load_total_kg")
blocks <- c("% change", trait_label(own_scale))
panel_data <- changes |>
  left_join(select(phenotype, subject, group), by = "subject") |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(block = if_else(trait %in% own_scale, trait_label(trait), blocks[1])) |>
  select(subject, group, trait, block, value)

pct_check <- with(
  phenotype, 100 * (leg_press_1rm_post_kg - leg_press_1rm_pre_kg) / leg_press_1rm_pre_kg
)
stopifnot(
  n_distinct(panel_data$subject) == 16,
  setequal(panel_data$group, c("HR", "LR")),
  all.equal(
    filter(panel_data, trait == "pct_leg_press_1rm")$value,
    pct_check[match(filter(panel_data, trait == "pct_leg_press_1rm")$subject, phenotype$subject)]
  ),
  all.equal(
    filter(panel_data, trait == "comp_hypertrophy")$value,
    phenotype$comp_hypertrophy[
      match(filter(panel_data, trait == "comp_hypertrophy")$subject, phenotype$subject)
    ]
  )
)

drawn <- panel_data |>
  mutate(
    label = factor(trait_label(trait, change = FALSE), rev(trait_label(trait_order))),
    block = factor(block, blocks),
    group = factor(group, c("HR", "LR"))
  )
stopifnot(!anyNA(drawn$label))
dodge <- position_dodge(width = 0.6, reverse = TRUE)
panel <- ggplot(drawn, aes(value, label, colour = group, shape = group)) +
  geom_vline(
    aes(xintercept = 0),
    data = distinct(filter(drawn, block != blocks[3]), block), linewidth = 0.2, colour = "grey60"
  ) +
  geom_point(position = dodge, size = 1, alpha = 0.9) +
  stat_summary(
    aes(group = group),
    fun = median, geom = "point", shape = "|", size = 3,
    position = dodge, show.legend = FALSE
  ) +
  ggh4x::facet_grid2(block ~ ., scales = "free", space = "free_y", independent = "x") +
  scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_colour_manual(values = group_greys, name = NULL) +
  scale_shape_manual(values = group_shapes, name = NULL) +
  labs(
    title = "Each person's training response",
    subtitle = "16 people, T1 to T2 (% change); bar: group median",
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(
    legend.position = "bottom", legend.margin = margin(t = -4),
    strip.text = element_blank(), panel.spacing.y = unit(2, "mm")
  )

ggsave(
  here("04_Figures", "F4", "b_reports", "panels", "a_trait_changes.pdf"), panel,
  width = 75, height = 70, units = "mm", device = cairo_pdf
)
