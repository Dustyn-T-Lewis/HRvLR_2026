library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))
source(here("04_Figures", "style", "module_labels.R"))

inputs <- c(characterisation = characterisation_input)
labels <- read_module_labels()
ora <- readxl::read_excel(characterisation_file, "ora")
top_members <- readxl::read_excel(characterisation_file, "top_members")
n_sets <- 3
holds <- paste(
  "Panel a: top three ORA sets per module (part ora) and the ten highest-kME members with their",
  "grid position (part member)."
)
n_columns <- 4

# GO Slim themes repeat GO BP names, so sets are counted by their printed label.
sets <- ora |>
  mutate(label = pathway_label(pathway, 40)) |>
  arrange(module, p) |>
  distinct(module, tolower(label), .keep_all = TRUE) |>
  slice_head(n = n_sets, by = module) |>
  mutate(row = row_number(), .by = module) |>
  select(module, row, set_id, collection, pathway, label, n_overlap, gene_ratio, p, fdr)
members <- top_members |>
  arrange(module, desc(kme)) |>
  mutate(rank = row_number(), .by = module) |>
  mutate(row = (rank - 1) %/% n_columns + 1, column = (rank - 1) %% n_columns + 1)

panel_data <- bind_rows(
  mutate(sets, part = "ora"),
  select(members, module, row, column, rank, uniprot_id, gene, kme) |> mutate(part = "member")
) |>
  relocate(part, module, row, column, rank) |>
  arrange(match(module, labels$module), part != "ora", row, column)

check <- inner_join(sets, ora, by = c("module", "set_id"))
stopifnot(
  setequal(sets$module, labels$module), all(count(sets, module)$n == n_sets),
  nrow(check) == nrow(sets), all(check$p.x == check$p.y), all(check$fdr.x == check$fdr.y),
  setequal(members$module, labels$module), max(members$row) <= n_sets
)

module_levels <- labels$module
facet_strips <- ggh4x::strip_themed(
  background_y = ggh4x::elem_list_rect(
    fill = scales::alpha(module_levels, 0.45), colour = "grey40", linewidth = 0.3
  ),
  text_y = ggh4x::elem_list_text(colour = "black", face = "bold")
)
fdr_cap <- 20
ora_plot <- sets |>
  mutate(
    module = factor(module, module_levels),
    row_id = factor(paste(module, row), rev(paste(rep(module_levels, each = n_sets), 1:n_sets))),
    log_fdr = pmin(-log10(fdr), fdr_cap)
  ) |>
  ggplot(aes(log_fdr, row_id)) +
  geom_vline(
    xintercept = -log10(alpha_fdr), linetype = "dashed", linewidth = 0.25,
    colour = "grey45"
  ) +
  geom_point(aes(fill = collection, size = n_overlap),
    shape = 21, colour = "grey20",
    stroke = 0.2
  ) +
  scale_fill_manual(
    values = pathway_collections, breaks = names(pathway_collections),
    labels = collection_names[names(pathway_collections)], name = "Database"
  ) +
  scale_size_area(max_size = 2.6, name = "members\nin set") +
  scale_y_discrete(labels = tibble::deframe(transmute(sets, paste(module, row), label))) +
  scale_x_continuous(
    limits = c(0, fdr_cap), breaks = c(0, 5, 10, 20), labels = c("0", "5", "10", "20+")
  ) +
  ggh4x::facet_grid2(
    module ~ .,
    scales = "free_y", space = "free_y", switch = "y", strip = facet_strips
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 2)), size = guide_legend(order = 2)
  ) +
  labs(
    title = "What each module holds",
    subtitle = sprintf(
      "Top %d ORA sets per module against the measured genes; dashed: FDR = %s",
      n_sets, alpha_fdr
    ),
    x = "-log10 FDR", y = NULL
  ) +
  figure_theme +
  theme(
    axis.text.y = element_text(size = 5), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 5.2, hjust = 0.5),
    panel.spacing.y = unit(0.6, "mm"), legend.position = "left",
    legend.text = element_text(size = 5.5), legend.title = element_text(size = 6)
  )

member_plot <- members |>
  mutate(module = factor(module, module_levels)) |>
  ggplot(aes(column, row)) +
  geom_tile(aes(fill = kme), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = gene), size = 1.8, fontface = "italic") +
  scale_fill_gradient(low = "grey96", high = "#9E9AC8", limits = c(NA, 1), name = "kME") +
  scale_y_reverse(expand = expansion(add = 0)) +
  scale_x_continuous(expand = expansion(add = 0)) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  guides(fill = guide_colourbar(barheight = unit(12, "mm"), barwidth = unit(2.5, "mm"))) +
  labs(title = " ", subtitle = "Ten highest-kME members", x = NULL, y = NULL) +
  figure_theme +
  theme(
    axis.text = element_blank(), panel.grid.major = element_blank(), strip.text = element_blank(),
    panel.spacing.y = unit(0.6, "mm"),
    legend.text = element_text(size = 5.5), legend.title = element_text(size = 6)
  )

panel <- (ora_plot | member_plot) + plot_layout(widths = c(1.15, 1))

ggsave(
  here("04_Figures", "S5", "b_reports", "panels", "a_characterisation.pdf"), panel,
  width = 180, height = 115, units = "mm", device = grDevices::cairo_pdf
)
