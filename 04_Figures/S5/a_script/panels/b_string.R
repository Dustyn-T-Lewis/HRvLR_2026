library(here)
library(dplyr)
library(ggplot2)
library(ggraph)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "module_labels.R"))

inputs <- c(characterisation = characterisation_input)
labels <- read_module_labels()
edges <- readxl::read_excel(characterisation_file, "network_edges")
string <- readxl::read_excel(characterisation_file, "string")
identity <- module_identity(labels, width = 20)
n_named <- 5
holds <- "Panel b: STRING edges drawn, among each module's 25 highest-kME members."

panel_data <- edges |>
  select(module, from_gene, to_gene, combined_score) |>
  arrange(match(module, labels$module), desc(combined_score))
nodes <- bind_rows(
  select(panel_data, module, gene = from_gene), select(panel_data, module, gene = to_gene)
) |>
  count(module, gene, name = "degree") |>
  mutate(named = rank(-degree, ties.method = "first") <= n_named, .by = module)

stopifnot(
  all(panel_data$module %in% labels$module), nrow(panel_data) == nrow(edges),
  all(panel_data$combined_score >= 700)
)

# Each module is laid out on its own, its components packed side by side (igraph
# layout_components), then all are drawn as one graph split by facet_nodes().
set.seed(20261015)
coordinates <- split(panel_data, panel_data$module) |>
  purrr::imap(\(e, m) {
    g <- igraph::graph_from_data_frame(select(e, from_gene, to_gene), directed = FALSE)
    xy <- igraph::layout_components(g)
    tibble(module = m, gene = igraph::V(g)$name, x = xy[, 1], y = xy[, 2])
  }) |>
  purrr::list_rbind()
graph <- tidygraph::tbl_graph(
  nodes = inner_join(nodes, coordinates, by = c("module", "gene")) |>
    mutate(name = paste(module, gene), module = factor(module, labels$module)),
  edges = transmute(
    panel_data,
    from = paste(module, from_gene), to = paste(module, to_gene), module, combined_score
  ),
  directed = FALSE, node_key = "name"
)

strip_text <- string |>
  filter(module %in% nodes$module) |>
  transmute(module, text = sprintf("%s\nobs/exp %.1f, %s", identity[module], ratio, label_fdr(fdr)))
panel <- ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(alpha = combined_score), colour = "grey45", edge_width = 0.2) +
  geom_node_point(aes(fill = I(as.character(module)), size = degree),
    shape = 21,
    colour = "grey25", stroke = 0.2
  ) +
  geom_node_text(aes(label = if_else(named, gene, NA_character_)),
    size = 1.8,
    fontface = "italic", repel = TRUE, max.overlaps = Inf, segment.size = 0.1,
    box.padding = 0.25, bg.colour = "white", bg.r = 0.1, na.rm = TRUE
  ) +
  facet_nodes(
    ~module,
    scales = "free", ncol = 5,
    labeller = as_labeller(tibble::deframe(strip_text))
  ) +
  scale_edge_alpha(range = c(0.25, 0.8), limits = c(700, 1000), name = "STRING score") +
  scale_size_area(max_size = 2.2, name = "edges") +
  labs(
    title = "STRING interactions within modules",
    subtitle = sprintf(
      paste(
        "Edges (combined score ≥ 700) among each module's 25 highest-kME members;",
        "the %d best-connected named. Strip: module-wide STRING enrichment"
      ),
      n_named
    )
  ) +
  figure_theme +
  theme(
    axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey85", linewidth = 0.3),
    strip.text = element_text(size = 5, lineheight = 0.95),
    legend.text = element_text(size = 5.5), legend.title = element_text(size = 6)
  )

ggsave(
  here("04_Figures", "S5", "b_reports", "panels", "b_string.pdf"), panel,
  width = 180, height = 95, units = "mm", device = grDevices::cairo_pdf
)
