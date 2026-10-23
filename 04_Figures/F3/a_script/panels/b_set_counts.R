# F3 b: gene sets per contrast at fgsea FDR < 0.05, up and down, over the four MSigDB databases
# (FDR within each). Read beside F1 d: the groups change a similar number of proteins, but HR moves
# many more sets.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))

inputs <- c(set_tests = "02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds")
holds <- "Panel b: per contrast and direction, gene sets at fgsea FDR < 0.05 and sets tested."
fgsea <- readRDS(here(inputs[["set_tests"]]))$set_tests |>
  filter(method == "fgsea", collection != "GO_Slim")
contrasts <- unique(fgsea$contrast)

# A set fgsea could not score (NES missing) counts as tested but in neither direction.
panel_data <- fgsea |>
  filter(fdr < alpha_fdr) |>
  count(contrast, direction = factor(if_else(nes > 0, "Up", "Down"), c("Up", "Down")),
    name = "n_fdr_05"
  ) |>
  complete(contrast = contrasts, direction, fill = list(n_fdr_05 = 0L)) |>
  inner_join(count(fgsea, contrast, name = "n_tested"), by = "contrast") |>
  mutate(
    family = contrast_panel(contrast), label = contrast_short(contrast),
    contrast = factor(contrast, rev(contrasts))
  ) |>
  select(family, label, contrast, direction, n_tested, n_fdr_05) |>
  arrange(family, desc(contrast), direction)

totals <- fgsea |>
  summarise(n = sum(fdr < alpha_fdr, na.rm = TRUE), .by = contrast) |>
  inner_join(summarise(panel_data, n_fdr_05 = sum(n_fdr_05), .by = contrast), by = "contrast")
stopifnot(
  nrow(totals) == length(contrasts), all(totals$n == totals$n_fdr_05),
  n_distinct(panel_data$n_tested) == 1
)

side <- \(direction) if_else(direction == "Up", 1, -1)
bars <- mutate(panel_data, signed = side(direction) * n_fdr_05)
limit <- max(panel_data$n_fdr_05) * 1.3
panel <- ggplot(bars, aes(signed, contrast)) +
  geom_col(aes(fill = direction), width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_text(
    aes(label = n_fdr_05, hjust = if_else(direction == "Up", -0.2, 1.2)),
    size = 2
  ) +
  scale_fill_manual(values = direction_colours, name = NULL) +
  scale_x_continuous(labels = abs, limits = c(-limit, limit)) +
  scale_y_discrete(labels = rlang::set_names(contrast_markup(contrasts), contrasts)) +
  ggh4x::facet_grid2(
    family ~ ., scales = "free_y", space = "free_y", switch = "y",
    strip = contrast_strips(y = levels(droplevels(bars$family))),
    labeller = labeller(family = label_wrap_gen(10))
  ) +
  labs(
    x = "gene sets (down | up)", y = NULL, title = "Gene sets per contrast",
    subtitle = sprintf(
      "fgsea FDR < %s within each database;\n%s sets tested (%s)",
      alpha_fdr, scales::comma(panel_data$n_tested[1]),
      paste(collection_names[setdiff(names(pathway_collections), "GO_Slim")], collapse = ", ")
    )
  ) +
  figure_theme +
  theme(
    axis.text.y.left = ggtext::element_markdown(size = 5.5), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 5.5), panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1, "mm"), legend.position = "bottom"
  )

save_panel(panel, "F3", "b_set_counts", 80, 110)
