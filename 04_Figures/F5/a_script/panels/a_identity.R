library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "module_labels.R"))

inputs <- c(
  characterisation = characterisation_input,
  modules = "01_Processing/06_Module_Building/c_data/modules.rds"
)
labels <- read_module_labels()
string <- readxl::read_excel(characterisation_file, "string")
modules <- readRDS(here(inputs[["modules"]]))
identity <- module_identity(labels)
holds <- paste(
  "Panel a: each module's size, top ORA set and STRING enrichment; identity is the label drawn."
)

panel_data <- labels |>
  select(module, n_proteins, top_set, top_set_fdr) |>
  mutate(identity = identity[module]) |>
  inner_join(
    select(string, module, n_edges, n_expected,
      string_ratio = ratio, string_p = p,
      string_fdr = fdr
    ),
    by = "module"
  )

stopifnot(
  setequal(panel_data$module, rownames(modules$eigengenes)),
  nrow(panel_data) == nrow(labels),
  all(panel_data$n_proteins == modules$module_summary$n_proteins[
    match(panel_data$module, modules$module_summary$module)
  ]),
  all(abs(panel_data$string_ratio - labels$string_ratio) < 1e-9)
)

d <- mutate(
  panel_data,
  module = factor(module, rev(module)),
  string_call = if_else(
    string_fdr < alpha_fdr, sprintf("FDR < %s", alpha_fdr), sprintf("FDR ≥ %s", alpha_fdr)
  )
)
size_plot <- ggplot(d, aes(n_proteins, module)) +
  geom_col(aes(fill = module), width = 0.65, colour = "grey30", linewidth = 0.15, alpha = 0.85) +
  geom_text(aes(label = n_proteins), hjust = -0.2, size = 1.8, colour = "grey30") +
  scale_fill_identity() +
  scale_y_discrete(labels = identity) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "Co-expression modules",
    subtitle = sprintf(
      "%d WGCNA modules, named by top ORA set (FDR < %s)\nSTRING dot filled at FDR < %s",
      nrow(d), alpha_fdr, alpha_fdr
    ),
    x = "proteins", y = NULL
  ) +
  figure_theme +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 5.5))
string_plot <- ggplot(d, aes(string_ratio, module)) +
  geom_vline(xintercept = 1, linewidth = 0.25, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 1, xend = string_ratio), linewidth = 0.3, colour = "grey60") +
  geom_point(aes(shape = string_call), size = 1.4, colour = "grey15", stroke = 0.4) +
  scale_shape_manual(
    values = rlang::set_names(c(16, 1), sprintf(c("FDR < %s", "FDR ≥ %s"), alpha_fdr)),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, NA), breaks = scales::breaks_width(2), expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = "STRING edges\nobs / exp", y = NULL) +
  figure_theme +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())

panel <- (size_plot | string_plot) + plot_layout(widths = c(3, 1.3))

ggsave(
  here("04_Figures", "F5", "b_reports", "panels", "a_identity.pdf"), panel,
  width = 100, height = 75, units = "mm", device = grDevices::cairo_pdf
)
