# What each module is: its enriched gene sets, its hub proteins, and whether its members
# interact more than chance in STRING.
#
# Enrichment is clusterProfiler's hypergeometric ORA against 03_Pathway_Enrichment's own sets,
# with the 1,900 measured genes as universe; testing against the genome would reward any module
# for being muscle. BH runs within each module.
#
# The STRING check is STRINGdb's PPI enrichment at combined score >= 700, with the background set
# to the measured proteins STRING can map. Its expected edge count accounts for each member's
# degree. Read from the local v12 files in 00_Input/downloads.
#
# The hub drawings show each module's 25 highest-kME members and the STRING edges among them.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(clusterProfiler)
  library(STRINGdb)
  library(patchwork)
})

stage <- here("04_Network", "02_characterise_modules")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
# Clear last run's figures so the bundle holds only this run's pages.
unlink(list.files(figure_dir, "[.](png|pdf)$", full.names = TRUE))

inputs <- c(
  modules = "04_Network/01_build_modules/c_data/modules.rds",
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  string_links = "00_Input/downloads/9606.protein.links.v12.0.txt.gz",
  string_aliases = "00_Input/downloads/9606.protein.aliases.v12.0.txt.gz",
  string_info = "00_Input/downloads/9606.protein.info.v12.0.txt.gz"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 01_build_modules and 00_build_gene_sets, and fetch STRING (00_Input/README.md). ",
    "Missing: ", paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
modules <- readRDS(paths[["modules"]])
gs <- readRDS(paths[["gene_sets"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)
members <- filter(modules$membership, module != "grey")
string_min_score <- 700L
hubs_per_module <- 10L

term2gene <- enframe(gs$sets, "set_id", "gene") |> unnest(gene)
ora_fit <- compareCluster(
  split(members$gene, members$module),
  fun = "enricher", TERM2GENE = term2gene, universe = gs$gene_universe,
  minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1
)
ora <- as_tibble(ora_fit) |>
  transmute(
    module = as.character(Cluster), set_id = ID, gene_ratio = GeneRatio, bg_ratio = BgRatio,
    count = Count, p = pvalue, fdr = p.adjust
  ) |>
  left_join(select(gs$set_catalog, set_id, database, pathway), by = "set_id")

hubs <- members |>
  slice_max(kme, n = hubs_per_module, by = module, with_ties = FALSE)

string_db <- STRINGdb$new(
  version = "12.0", species = 9606, score_threshold = string_min_score,
  network_type = "full", input_directory = here("00_Input", "downloads")
)
mapped <- string_db$map(
  as.data.frame(modules$membership), "uniprot_id",
  removeUnmappedRows = TRUE
)
string_db$set_background(mapped$STRING_id)
string_check <- mapped |>
  filter(module != "grey") |>
  summarise(n_mapped = n(), ids = list(STRING_id), .by = module) |>
  mutate(
    enrichment = map(ids, string_db$get_ppi_enrichment),
    edges = map_dbl(enrichment, "edges"),
    expected = map_dbl(enrichment, "lambda"),
    ratio = edges / expected,
    p = map_dbl(enrichment, "enrichment"),
    fdr = p.adjust(p, "BH")
  ) |>
  select(-ids, -enrichment)

labels <- count(members, module, name = "proteins") |>
  left_join(
    ora |>
      slice_min(p, n = 1, by = module, with_ties = FALSE) |>
      select(module, top_set = pathway, top_set_fdr = fdr),
    by = "module"
  ) |>
  left_join(
    ora |>
      filter(database == "GO_Slim") |>
      slice_min(p, n = 1, by = module, with_ties = FALSE) |>
      select(module, top_go_slim = pathway, go_slim_fdr = fdr),
    by = "module"
  ) |>
  left_join(summarise(hubs, hubs = paste(gene, collapse = ", "), .by = module), by = "module") |>
  left_join(select(string_check, module, string_ratio = ratio, string_fdr = fdr), by = "module") |>
  arrange(desc(proteins))
print(as.data.frame(select(labels, module, proteins, top_set, top_set_fdr, string_ratio)),
  digits = 3
)


# ---- figures -------------------------------------------------------------------------------

save_figure <- function(figure, name, width, height) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
      width = width, height = height, dpi = 200, bg = "white"
    )
  })
}
caption_theme <- theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0))

save_figure(
  ora_fit |>
    clusterProfiler::filter(p.adjust < 0.05) |>
    enrichplot::dotplot(showCategory = 3, label_format = 50, font.size = 7) +
    labs(
      title = "Module enrichment",
      subtitle = sprintf(
        "clusterProfiler ORA, universe %d measured genes; top 3 sets per module at BH < 0.05",
        length(gs$gene_universe)
      ),
      caption = paste(
        "enrichplot::dotplot of the compareCluster result. Columns: modules, members carrying a",
        "set in brackets. Size is gene ratio, colour BH within module; modules with no set at",
        "BH < 0.05 are absent. Table: c_data/02_characterise_modules.xlsx, ora."
      )
    ) +
    caption_theme,
  "01_module_enrichment", 10, 8
)

module_colours <- set_names(string_check$module, string_check$module)
save_figure(
  string_check |>
    mutate(module = factor(module, levels = rev(labels$module))) |>
    ggplot(aes(y = module)) +
    geom_segment(aes(x = expected, xend = edges, yend = module), colour = "grey70") +
    geom_point(aes(x = expected), shape = 124, size = 4) +
    geom_point(aes(x = edges, fill = module), shape = 21, size = 3) +
    geom_text(aes(x = edges, label = sprintf("%.1fx", ratio)), hjust = -0.4, size = 2.8) +
    scale_fill_manual(values = module_colours, guide = "none") +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      x = "within-module edges (log scale)", y = NULL,
      title = "STRING edges within each module",
      subtitle = sprintf(
        "STRINGdb PPI enrichment, v12, score >= %d, background the mapped measured proteins",
        string_min_score
      ),
      caption = paste(
        "Point: observed high-confidence edges among a module's members. Tick: edges STRINGdb",
        "expects from the members' degrees in the background. Label: observed over expected.",
        "Table: c_data/02_characterise_modules.xlsx, string."
      )
    ) +
    theme_minimal(base_size = 10) +
    caption_theme,
  "02_string_edges", 8, 5
)

# Each module's highest-kME members and the STRING edges among them, at the same score threshold
# as the enrichment test. Node fill is kME, size the node's degree among the drawn hubs. The
# layout is seeded; its geometry carries no meaning.
hubs_drawn <- 25L
hub_nodes <- mapped |>
  filter(module != "grey") |>
  slice_max(kme, n = hubs_drawn, by = module, with_ties = FALSE)
hub_edges <- string_db$get_interactions(hub_nodes$STRING_id) |>
  distinct(from, to, .keep_all = TRUE) |>
  inner_join(select(hub_nodes, from = STRING_id, module), by = "from") |>
  inner_join(select(hub_nodes, to = STRING_id, to_module = module), by = "to") |>
  filter(module == to_module) |>
  select(module, from, to, combined_score)
max_degree <- max(table(c(hub_edges$from, hub_edges$to)))
draw_hub_network <- function(which_module) {
  nodes <- hub_nodes |>
    filter(module == which_module) |>
    transmute(name = STRING_id, gene, kme)
  edges <- hub_edges |>
    filter(module == which_module) |>
    transmute(from, to, score = combined_score / 1000)
  graph <- nodes |>
    tidygraph::tbl_graph(edges = edges, directed = FALSE, node_key = "name") |>
    tidygraph::activate(nodes) |>
    mutate(degree = tidygraph::centrality_degree()) |>
    filter(degree > 0)
  isolated <- nrow(nodes) - igraph::vcount(graph)
  set.seed(42)
  ggraph::ggraph(graph, layout = "fr") +
    ggraph::geom_edge_link(aes(edge_width = score), colour = "grey70", alpha = 0.6) +
    ggraph::geom_node_point(aes(size = degree, fill = kme), shape = 21, colour = "grey30") +
    ggraph::geom_node_text(aes(label = gene), size = 2.2, repel = TRUE, max.overlaps = Inf) +
    ggraph::scale_edge_width(range = c(0.2, 1.2), guide = "none") +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0.4, 1), name = "kME") +
    scale_size_continuous(range = c(1.5, 5), limits = c(1, max_degree), name = "degree") +
    labs(title = sprintf(
      "%s: %d edges among %d hubs, %d without an edge", which_module, nrow(edges), nrow(nodes),
      isolated
    )) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10), legend.position = "right")
}
hub_pages <- labels$module |>
  map(draw_hub_network) |>
  split(ceiling(seq_along(labels$module) / 4)) |>
  imap(\(panels, page) {
    wrap_plots(panels, ncol = 2, guides = "collect") +
      plot_annotation(
        title = "Module hub networks",
        subtitle = sprintf(
          "top %d members by kME per module, STRING v12 edges at score >= %d; page %s of %d",
          hubs_drawn, string_min_score, page, ceiling(length(labels$module) / 4)
        ),
        caption = stringr::str_wrap(width = 190, paste(
          "Nodes: a module's highest-kME members with at least one edge among them, fill = kME,",
          "size = degree among the drawn hubs; the title counts those without an edge. Edges:",
          "STRING combined score, width = score. Layout is force-directed and seeded; its",
          "geometry carries no meaning. Table: c_data/02_characterise_modules.xlsx, hub_edges."
        )),
        theme = theme(
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey40")
        )
      )
  })
# Paged, so PDF only.
pdf(file.path(figure_dir, "03_hub_networks.pdf"), width = 11, height = 10, bg = "white")
walk(hub_pages, print)
invisible(dev.off())
message("drew 03_hub_networks: ", length(hub_pages), " pages")

packages <- c(
  "here", "clusterProfiler", "enrichplot", "STRINGdb", "tidygraph", "ggraph", "dplyr", "purrr",
  "ggplot2"
)
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    labels = labels, ora = ora, ora_fit = ora_fit, hubs = hubs, string = string_check,
    hub_edges = hub_edges,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "module_characterisation.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    labels = labels, ora = arrange(ora, module, p), hubs = hubs, string = string_check,
    hub_edges = left_join(hub_edges,
      select(hub_nodes, from = STRING_id, from_gene = gene),
      by = "from"
    ) |>
      left_join(select(hub_nodes, to = STRING_id, to_gene = gene), by = "to"),
    input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "02_characterise_modules.xlsx")
)
combined <- file.path(figure_dir, "02_characterise_modules_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote module_characterisation.rds, 02_characterise_modules.xlsx and ",
  length(pages), " figures bundled into ", qpdf::pdf_length(combined), " pages"
)
