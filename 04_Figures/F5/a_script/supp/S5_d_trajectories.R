library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "module_labels.R"))

inputs <- c(
  modules = "01_Processing/06_Module_Building/c_data/module_building.rds",
  module_tests = "02_HRvLR/03_Module_Contrasts/c_data/module_contrasts.rds",
  characterisation = characterisation_input
)
modules <- readRDS(here(inputs[["modules"]]))
tests <- readRDS(here(inputs[["module_tests"]]))$eigengene_tests
identity <- module_identity(read_module_labels())
groups <- c("HR", "LR")

# The module drawn is the one with the smallest eigengene p in a training, acute or HR - LR
# contrast; the interactions are left out because they test a difference of changes, which the
# per-person lines show less directly.
eligible <- filter(tests, contrast_panel(contrast) != "Interaction")
chosen <- slice_min(eligible, p, n = 1, with_ties = FALSE)
module <- chosen$module

panel_data <- tibble(
  sample_id = colnames(modules$eigengenes), eigengene = modules$eigengenes[module, ]
) |>
  inner_join(select(modules$meta, sample_id, subject, group, timepoint), by = "sample_id") |>
  mutate(module = module, chosen_contrast = chosen$contrast, chosen_p = chosen$p) |>
  arrange(group, subject, timepoint)

stopifnot(
  nrow(panel_data) == ncol(modules$eigengenes), !anyNA(panel_data$eigengene),
  setequal(panel_data$group, groups),
  all(table(panel_data$subject) == n_distinct(panel_data$timepoint)), chosen$p == min(eligible$p)
)

d <- mutate(panel_data, group = factor(group, groups))
means <- summarise(d, eigengene = mean(eigengene), .by = c(group, timepoint))
panel <- ggplot(d, aes(timepoint, eigengene, colour = group)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_line(aes(group = subject), alpha = 0.6, linewidth = 0.25) +
  geom_line(aes(group = group), data = means, linewidth = 1) +
  geom_point(data = means, size = 1.4) +
  facet_wrap(~group, nrow = 1) +
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.2)) +
  labs(
    title = sprintf("%s, person by person", identity[[module]]),
    subtitle = sprintf(
      "Exploratory: the smallest eigengene p outside the interactions\n%s logFC %.2f, %s, %s\n%s",
      contrast_label(chosen$contrast), chosen$logFC, label_p(chosen$p), label_fdr(chosen$fdr),
      "Thin lines: people; thick: group mean"
    ),
    x = NULL, y = "eigengene"
  ) +
  figure_theme +
  theme(strip.text = element_text(face = "bold", size = 6.5))

save_panel(panel, "S5", "d_trajectories", 110, 95)
