# F5 b: after the extra bout, which modules change together across people, within HR and within
# LR (02_HRvLR/08_Module_Detail). Boxed: pairs whose HR and LR correlations differ at Fisher
# p < 0.05. The training window is in 04_Figures/S7.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "module_labels.R"))

inputs <- c(
  module_detail = "02_HRvLR/08_Module_Detail/c_data/module_detail.rds",
  characterisation = characterisation_input
)
labels <- read_module_labels()
names_short <- unname(module_identity(labels, width = 26)[labels$module])
panel_data <- readRDS(here(inputs[["module_detail"]]))$coupling |>
  filter(window == "acute") |>
  select(module_a, module_b, rho_HR, rho_LR, difference, p_difference)
stopifnot(nrow(panel_data) == choose(nrow(labels), 2))

d <- panel_data |>
  pivot_longer(c(rho_HR, rho_LR), names_to = "group", values_to = "rho") |>
  mutate(
    group = factor(group, c("rho_HR", "rho_LR"), c("within HR", "within LR")),
    first = match(module_a, labels$module), second = match(module_b, labels$module),
    x = factor(names_short[pmin(first, second)], names_short),
    y = factor(names_short[pmax(first, second)], rev(names_short)),
    differs = p_difference < alpha_nominal
  )
panel <- ggplot(d, aes(x, y, fill = rho)) +
  geom_tile(colour = "white") +
  geom_tile(data = filter(d, differs), fill = NA, colour = "black", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", rho)), size = 1.7) +
  scale_fill_gradient2(
    low = "#4393C3", high = "#D6604D", limits = c(-1, 1), name = "Spearman\nrho"
  ) +
  facet_grid(. ~ group) +
  labs(
    title = "After the bout, HR's modules move together",
    subtitle = sprintf(
      paste(
        "Spearman rho between two modules' acute eigengene changes across 8 people;",
        "box: HR and LR differ at Fisher p < %s (%d of %d pairs)"
      ),
      alpha_nominal, sum(panel_data$p_difference < alpha_nominal), nrow(panel_data)
    ),
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank(), strip.text = element_text(face = "bold")
  )
save_panel(panel, "F5", "b_coupling", 140, 85)
