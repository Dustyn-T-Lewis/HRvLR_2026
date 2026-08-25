# What are the twelve modules?
#
# A module hit that cannot be named is not a result. This assigns each module a
# biological description from GO over-representation across all three
# sub-ontologies, and lists the proteins that define it by module membership.
#
# The universe is the measured proteome, not the human genome. That single
# choice is what separates this from the two enrichment nulls this project has
# already found miscalibrated: fgsea's preranked null permutes gene labels and
# so treats co-regulated proteins as exchangeable, and STRING's PPI p was
# invalid on an MS proteome. A hypergeometric test against the proteins that
# were actually detected asks a well-posed question instead.
#
# This is labelling, not evidence. Nothing downstream gates on an enrichment
# q-value, and V1 recorded an ORA hit computed against a correct universe that
# later failed leave-one-subject-out.

pacman::p_load(
  here, dplyr, purrr, tibble, readr, clusterProfiler, org.Hs.eg.db, openxlsx
)

source(here("functions", "association.R"))

OUT_DIR <- here("03_Features", "c_data")
TOP_HUBS <- 15

membership <- openxlsx::read.xlsx(
  file.path(OUT_DIR, "01_modules.xlsx"), "module_membership"
) |>
  filter(.data$module != "grey", !is.na(.data$gene), .data$gene != "")

universe <- unique(membership$gene)

# Module membership: how strongly each protein correlates with its own module's
# eigengene. The top of that ranking is what a module is made of, and it does
# not depend on any enrichment database being right.
imputed <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
abund <- as.matrix(imputed$data)
rownames(abund) <- imputed$annotation$gene

eigengenes <- module_matrix()

kme <- map_dfr(rownames(eigengenes), function(me) {
  mod <- sub("^ME_", "", me)
  genes <- membership$gene[membership$module == mod]
  genes <- intersect(genes, rownames(abund))
  tibble(
    module = mod,
    gene = genes,
    kme = as.numeric(stats::cor(
      t(abund[genes, colnames(eigengenes), drop = FALSE]),
      eigengenes[me, ]
    ))
  )
}) |>
  arrange(.data$module, desc(.data$kme))

hubs <- kme |>
  slice_max(.data$kme, n = TOP_HUBS, by = "module") |>
  summarise(
    hub_proteins = paste(.data$gene, collapse = ", "),
    .by = "module"
  )

ontology_hits <- function(genes, ont) {
  res <- clusterProfiler::enrichGO(
    gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL", ont = ont, pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = FALSE
  )
  if (is.null(res) || !nrow(as.data.frame(res))) {
    return(tibble())
  }
  as_tibble(as.data.frame(res)) |>
    mutate(ontology = ont, .before = 1)
}

set.seed(42)
enrichment <- map_dfr(sort(unique(membership$module)), function(mod) {
  genes <- membership$gene[membership$module == mod]
  map_dfr(c("BP", "CC", "MF"), function(ont) {
    ontology_hits(genes, ont) |> mutate(module = mod, .before = 1)
  })
})

# One line per module: its size, its top term in each sub-ontology, and its
# hubs. This is the table a reader consults when a module name appears in a
# result.
top_term <- function(mod, ont) {
  d <- enrichment |>
    filter(.data$module == mod, .data$ontology == ont) |>
    slice_min(.data$p.adjust, n = 1, with_ties = FALSE)
  if (!nrow(d)) {
    return(NA_character_)
  }
  sprintf("%s (q=%.1e)", d$Description, d$p.adjust)
}

atlas <- tibble(module = sort(unique(membership$module))) |>
  mutate(
    n_proteins = map_int(.data$module, ~ sum(membership$module == .x)),
    n_terms = map_int(.data$module, ~ sum(enrichment$module == .x)),
    BP = map_chr(.data$module, top_term, ont = "BP"),
    CC = map_chr(.data$module, top_term, ont = "CC"),
    MF = map_chr(.data$module, top_term, ont = "MF")
  ) |>
  left_join(hubs, by = "module") |>
  arrange(desc(.data$n_proteins))

write.xlsx(
  list(
    module_atlas = atlas, enrichment = enrichment, module_membership_kme = kme
  ),
  file.path(OUT_DIR, "04_module_annotation.xlsx")
)

print(as.data.frame(atlas[, c("module", "n_proteins", "n_terms", "BP")]),
  right = FALSE
)
message(
  "\n", nrow(enrichment), " enriched terms across ",
  dplyr::n_distinct(enrichment$module), " of ", nrow(atlas),
  " modules, universe = ", length(universe), " detected proteins"
)
