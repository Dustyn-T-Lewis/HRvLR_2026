library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "module_labels.R"))

inputs <- c(
  module_tests = "02_HRvLR/03_Module_Contrasts/c_data/module_tests.rds",
  characterisation = characterisation_input
)
tests <- readRDS(here(inputs[["module_tests"]]))
labels <- read_module_labels()
contrasts <- unique(tests$eigengene_tests$contrast)
holds <- paste(
  "Panel b: eigengene limma results and fry on the module's proteins, every module and",
  "contrast."
)

panel_data <- tests$eigengene_tests |>
  select(contrast, module, logFC, logFC_lower, logFC_upper, t, p, fdr) |>
  left_join(
    select(tests$fry_tests, contrast, module,
      fry_direction = direction, fry_p = p,
      fry_fdr = fdr
    ),
    by = c("contrast", "module")
  ) |>
  arrange(match(module, labels$module), match(contrast, contrasts))

check <- inner_join(panel_data, tests$fry_tests, by = c("contrast", "module"))
stopifnot(
  length(contrasts) == 9, setequal(panel_data$module, labels$module),
  nrow(panel_data) == 9 * nrow(labels), nrow(check) == nrow(panel_data),
  all(check$fry_p == check$p.y), !anyNA(panel_data$logFC)
)

p_cap <- 3
calls <- sprintf(c("p < %s", "p ≥ %s"), alpha_nominal)
marks <- c(sprintf("FDR < %s", alpha_fdr), sprintf("fry p < %s", alpha_nominal))
d <- panel_data |>
  mutate(
    family = contrast_panel(contrast),
    column = factor(contrast, contrasts),
    module = factor(module, rev(labels$module)),
    log_p = pmin(-log10(p), p_cap),
    call = factor(if_else(p < alpha_nominal, calls[1], calls[2]), calls)
  )
limit <- ceiling(max(abs(d$logFC)) * 20) / 20

panel <- ggplot(d, aes(column, module)) +
  geom_point(aes(size = log_p, colour = logFC, shape = call), stroke = 0.45) +
  geom_point(
    aes(shape = mark), filter(d, fry_p < alpha_nominal) |> mutate(mark = marks[2]),
    size = 3.4, colour = "grey45", stroke = 0.35
  ) +
  geom_point(
    aes(size = log_p, shape = mark), filter(d, fdr < alpha_fdr) |> mutate(mark = marks[1]),
    colour = "black", stroke = 0.7
  ) +
  scale_colour_distiller(
    palette = "PiYG", direction = 1, limits = c(-limit, limit), name = "eigengene\nlogFC"
  ) +
  scale_size_area(
    max_size = 2.6, limits = c(0, p_cap), breaks = c(1, 2, p_cap),
    labels = c("1", "2", paste0(p_cap, "+")), name = "-log10 p"
  ) +
  scale_shape_manual(
    values = rlang::set_names(c(16, 1, 1, 0), c(calls, marks)), limits = c(calls, marks),
    name = NULL, drop = FALSE
  ) +
  scale_x_discrete(labels = rlang::set_names(contrast_markup(contrasts), contrasts)) +
  guides(
    colour = guide_colourbar(order = 1, barheight = unit(14, "mm"), barwidth = unit(2.5, "mm")),
    size = guide_legend(order = 2),
    shape = guide_legend(order = 3, override.aes = list(
      size = c(1.6, 1.6, 1.6, 3), colour = c("grey30", "grey30", "black", "grey45"),
      stroke = c(0.45, 0.45, 0.7, 0.35)
    ))
  ) +
  ggh4x::facet_grid2(
    . ~ family,
    scales = "free_x", space = "free_x", strip = contrast_strips(x = levels(d$family)),
    labeller = labeller(family = label_wrap_gen(10))
  ) +
  labs(
    x = NULL, y = NULL, title = "Module eigengenes across the nine contrasts",
    subtitle = sprintf(
      paste0(
        "limma on each eigengene; square: fry on the module's proteins\n",
        "p < %s in %d of %d tests, %.1f by chance"
      ),
      alpha_nominal, sum(d$p < alpha_nominal), nrow(d), alpha_nominal * nrow(d)
    )
  ) +
  figure_theme +
  matrix_theme +
  theme(axis.text.y = element_text(size = 5.5))

ggsave(
  here("04_Figures", "F5", "b_reports", "panels", "b_eigengene_matrix.pdf"), panel,
  width = 80, height = 75, units = "mm", device = grDevices::cairo_pdf
)
