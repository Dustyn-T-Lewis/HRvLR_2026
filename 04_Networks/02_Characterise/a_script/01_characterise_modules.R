# What each module is: its enriched gene sets, its hub proteins, its GO-Slim
# theme, and whether its members interact more than chance in STRING.
#
# Enrichment is clusterProfiler's hypergeometric ORA against the pathway
# stage's own gene sets, with the 1900 detected proteins as universe; testing
# against the genome would reward any module for being muscle. BH runs within
# each module.
#
# The STRING check asks whether a module's members share more high-confidence
# edges (combined score >= 700) than the same number of proteins drawn from
# this proteome at random. The null shuffles module labels over the proteins
# STRING can map, so it inherits this dataset's detection bias rather than
# STRING's own background.

pacman::p_load(here, dplyr, tidyr, tibble, readr, purrr, clusterProfiler)

source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_utils.R"))

set.seed(42)

OUT_DIR <- here("04_Networks", "02_Characterise", "c_data")
STRING_MIN_SCORE <- 700L
N_PERM <- 999L
N_HUBS <- 10L

modules <- readRDS(here("04_Networks", "01_Modules", "c_data", "modules.rds"))
gs <- readRDS(here("03_Pathways", "01_Gene_Sets", "c_data", "gene_sets.rds"))
members <- filter(modules$membership, .data$module != "grey")

term2gene <- enframe(gs$sets, "set", "uniprot_id") |> unnest("uniprot_id")
ora <- compareCluster(
  split(members$uniprot_id, members$module),
  fun = "enricher", TERM2GENE = term2gene, universe = gs$universe,
  minGSSize = SET_FLOOR, maxGSSize = SET_CEILING,
  pvalueCutoff = 1, qvalueCutoff = 1
) |>
  as_tibble() |>
  transmute(
    module = as.character(.data$Cluster), set = .data$ID,
    gene_ratio = .data$GeneRatio, bg_ratio = .data$BgRatio,
    count = .data$Count, p = .data$pvalue, bh = .data$p.adjust
  ) |>
  left_join(select(gs$catalog, "set", "collection", "theme"), by = "set")

hubs <- members |>
  slice_max(.data$kme, n = N_HUBS, by = "module", with_ties = FALSE)

module_theme <- ora |>
  filter(.data$bh < 0.05, !is.na(.data$theme)) |>
  count(.data$module, .data$theme, name = "n_sets") |>
  slice_max(.data$n_sets, n = 1, by = "module", with_ties = FALSE)

aliases <- read_tsv(
  here("00_input", "downloads", "9606.protein.aliases.v12.0.txt.gz"),
  col_names = c("string_id", "alias", "source"), comment = "#",
  show_col_types = FALSE
) |>
  filter(
    .data$source == "UniProt_AC",
    .data$alias %in% modules$membership$uniprot_id
  ) |>
  distinct(.data$alias, .keep_all = TRUE)
links <- read_delim(
  here("00_input", "downloads", "9606.protein.links.v12.0.txt.gz"),
  delim = " ", show_col_types = FALSE
) |>
  filter(
    .data$combined_score >= STRING_MIN_SCORE,
    .data$protein1 %in% aliases$string_id,
    .data$protein2 %in% aliases$string_id,
    .data$protein1 < .data$protein2
  )

mapped <- modules$membership |>
  inner_join(aliases, by = c(uniprot_id = "alias"))
label <- setNames(mapped$module, mapped$string_id)
edges_i <- links$protein1
edges_j <- links$protein2

within_edges <- function(lab) {
  same <- lab[edges_i] == lab[edges_j]
  table(factor(lab[edges_i][same], levels = unique(members$module)))
}
observed <- within_edges(label)
null <- replicate(N_PERM, within_edges(setNames(sample(label), names(label))))

string_check <- tibble(
  module = names(observed),
  n_mapped = as.integer(table(label)[names(observed)]),
  edges = as.integer(observed),
  null_median = apply(null, 1, stats::median),
  null_lo = apply(null, 1, stats::quantile, 0.025),
  null_hi = apply(null, 1, stats::quantile, 0.975),
  enrichment = .data$edges / .data$null_median,
  p = (1 + rowSums(null >= as.integer(observed))) / (N_PERM + 1)
) |>
  mutate(bh = p.adjust(.data$p, "BH"))

top_set <- ora |>
  slice_min(.data$p, n = 1, by = "module", with_ties = FALSE) |>
  select("module", top_set = "set", top_set_bh = "bh")
labels <- count(members, .data$module, name = "n_proteins") |>
  left_join(top_set, by = "module") |>
  left_join(select(module_theme, "module", "theme"), by = "module") |>
  left_join(
    summarise(hubs, hubs = paste(.data$gene, collapse = ", "), .by = "module"),
    by = "module"
  ) |>
  left_join(
    select(string_check, "module",
      string_enrichment = "enrichment", string_p = "p"
    ),
    by = "module"
  )

clear_dir(OUT_DIR)
characterisation <- list(
  labels = labels, ora = ora, hubs = hubs, string = string_check,
  string_min_score = STRING_MIN_SCORE, n_perm = N_PERM
)
saveRDS(characterisation, file.path(OUT_DIR, "module_characterisation.rds"))
write_workbook(
  file.path(OUT_DIR, "module_characterisation.xlsx"),
  characterisation[c("labels", "ora", "hubs", "string")]
)

print(as.data.frame(select(labels, -"hubs")), row.names = FALSE, digits = 3)
