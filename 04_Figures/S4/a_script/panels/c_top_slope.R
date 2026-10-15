# The strongest slope difference (smallest p across every feature type, window and trait), person
# by person: the feature's value in the window against the trait, a least-squares line per group,
# and each group's Spearman rho from 02_HRvLR/07_Group_Trait_Slopes.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("02_HRvLR", "helpers.R"))

inputs <- c(
  group_trait_slopes = "02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds",
  corrected = "01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds",
  phenotype = "00_Setup/Input/phenotype.csv", helpers = "02_HRvLR/helpers.R"
)
slopes <- readRDS(here(inputs[["group_trait_slopes"]]))
top <- slice_min(slopes$group_difference, p, n = 1, with_ties = FALSE)
holds <- "Panel c: the smallest-p feature's value and trait per person."

feature_files <- c(
  sets = "01_Processing/05_Set_Scores/c_data/singscore.rds",
  modules = "01_Processing/06_Module_Building/c_data/modules.rds"
)
if (top$feature_type != "proteins") {
  inputs <- c(inputs, setNames(feature_files[[top$feature_type]], top$feature_type))
}
corrected <- readRDS(here(inputs[["corrected"]]))
meta <- as_tibble(corrected$targets)
values <- switch(top$feature_type,
  proteins = corrected$E,
  sets = readRDS(here(inputs[["sets"]]))$scores,
  modules = readRDS(here(inputs[["modules"]]))$eigengenes
)
stopifnot(identical(colnames(values), meta$sample_id))
phenotype <- readr::read_csv(here(inputs[["phenotype"]]), show_col_types = FALSE)
traits <- trait_windows(phenotype)[[top$window]]
person <- window_values(values[top$feature, , drop = FALSE], meta)[[top$window]][1, ]

panel_data <- tibble(subject = names(person), feature_value = unname(person)) |>
  inner_join(select(traits, subject, trait_value = all_of(top$trait)), by = "subject") |>
  left_join(distinct(meta, subject, group), by = "subject") |>
  filter(!is.na(trait_value)) |>
  mutate(
    feature_type = top$feature_type, feature = top$feature, label = top$label,
    window = top$window, trait = top$trait, .before = 1
  )

rho <- slopes$within_group |>
  filter(
    feature_type == top$feature_type, feature == top$feature, window == top$window,
    trait == top$trait
  )
rho_check <- panel_data |>
  summarise(rho = cor(feature_value, trait_value, method = "spearman"), n = n(), .by = group)
stopifnot(
  nrow(panel_data) == top$n_people,
  isTRUE(all.equal(unname(rho$rho), rho_check$rho[match(rho$group, rho_check$group)])),
  all(rho$n_people == rho_check$n[match(rho$group, rho_check$group)])
)

value_name <- c(proteins = "log2 abundance", sets = "set score", modules = "eigengene")
feature_name <- if (top$feature_type == "sets") pathway_label(top$label, 45) else top$label
legend_text <- with(rho, setNames(sprintf("%s: rho %.2f (n = %d)", group, rho, n_people), group))
drawn <- mutate(panel_data, group = factor(group, c("HR", "LR")))

panel <- ggplot(drawn, aes(feature_value, trait_value, colour = group, shape = group)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = group_greys, labels = legend_text, name = NULL) +
  scale_shape_manual(values = group_shapes, labels = legend_text, name = NULL) +
  labs(
    title = sprintf("%s, person by person", feature_name),
    subtitle = sprintf(
      "HR slope minus LR slope: t = %.1f, %s, %s", top$t, label_p(top$p), label_fdr(top$fdr)
    ),
    x = sprintf("%s, %s", value_name[[top$feature_type]], screen_windows[[top$window]]),
    y = trait_label(top$trait)
  ) +
  figure_theme +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.margin = margin(t = -4))

ggsave(
  here("04_Figures", "S4", "b_reports", "panels", "c_top_slope.pdf"), panel,
  width = 85, height = 90, units = "mm", device = cairo_pdf
)
