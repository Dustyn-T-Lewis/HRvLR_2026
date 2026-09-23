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

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(clusterProfiler)
  library(STRINGdb)
})

stage <- here("04_Network", "02_characterise_modules")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

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

packages <- c("here", "clusterProfiler", "enrichplot", "STRINGdb", "dplyr", "purrr", "ggplot2")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    labels = labels, ora = ora, ora_fit = ora_fit, hubs = hubs, string = string_check,
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
    input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "02_characterise_modules.xlsx")
)
combined <- file.path(figure_dir, "02_characterise_modules_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote module_characterisation.rds, 02_characterise_modules.xlsx and a ",
  length(pages), "-page figure PDF"
)
