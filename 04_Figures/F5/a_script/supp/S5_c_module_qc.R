library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(
  modules = "01_Processing/06_Module_Building/c_data/module_building.rds",
  module_building = "01_Processing/06_Module_Building/c_data/module_building.xlsx"
)
modules <- readRDS(here(inputs[["modules"]]))
module_qc <- read_sheet(inputs[["module_building"]], "module_qc")
soft_threshold <- read_sheet(inputs[["module_building"]], "soft_threshold")
chosen_power <- soft_threshold$Power[soft_threshold$chosen]
other_powers <- sort(unique(as.numeric(module_qc$power_other)))

networks <- module_qc |>
  distinct(power = as.numeric(power_other), n_modules = n_modules_other, n_grey = grey_other) |>
  add_row(
    power = chosen_power, n_modules = nrow(modules$eigengenes),
    n_grey = sum(modules$membership$module == "grey")
  ) |>
  arrange(power)
kept <- module_qc |>
  transmute(power = as.numeric(power_other), module, n_proteins, kept_share)
panel_data <- bind_rows(
  mutate(networks, part = "network"), mutate(kept, part = "kept_share")
) |>
  relocate(part)

stopifnot(
  length(chosen_power) == 1, setequal(networks$power, c(other_powers, chosen_power)),
  anyDuplicated(networks$power) == 0,
  setequal(kept$module, rownames(modules$eigengenes)),
  all(kept$kept_share >= 0 & kept$kept_share <= 1)
)

power_axis <- scale_x_continuous(breaks = networks$power)
chosen_line <- geom_vline(
  xintercept = chosen_power, linewidth = 0.3, linetype = "dashed", colour = "grey45"
)
count_plot <- function(column, title, y) {
  ggplot(networks, aes(power, .data[[column]])) +
    chosen_line +
    geom_line(linewidth = 0.3, colour = "grey50") +
    geom_point(aes(colour = power == chosen_power), size = 1.4) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey55"), guide = "none") +
    power_axis +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, x = "soft power", y = y) +
    figure_theme
}
modules_plot <- count_plot("n_modules", "Modules per soft power", "modules") +
  labs(subtitle = sprintf("Power %d used (dashed)", chosen_power))
grey_plot <- count_plot("n_grey", "Proteins in no module", "grey proteins") +
  labs(subtitle = sprintf("of %s measured", scales::comma(nrow(modules$membership))))
kept_plot <- ggplot(kept, aes(power, kept_share)) +
  geom_line(aes(group = module, colour = module), linewidth = 0.3, alpha = 0.7) +
  geom_point(aes(fill = module), shape = 21, colour = "grey25", stroke = 0.2, size = 1.3) +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_x_continuous(breaks = sort(unique(kept$power))) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(
    title = sprintf("Power-%d modules at other powers", chosen_power),
    subtitle = "Share of each module's proteins in its best match",
    x = "soft power", y = "kept"
  ) +
  figure_theme

panel <- modules_plot / grey_plot / kept_plot

save_panel(panel, "S5", "c_module_qc", 95, 115)
