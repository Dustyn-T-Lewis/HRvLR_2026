# F04 Panel B: pathway over-representation by WGCNA module.
# Universe = full measured proteome (every protein in the network); bars filled by
# module colour, names inside, free x per facet, source database on the border.
pacman::p_load(here, dplyr, ggplot2, stringr, shadowtext, clusterProfiler, org.Hs.eg.db)
if (!exists("wgcna_mem")) source(here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))

wgcna_non_grey <- wgcna_mem |> filter(group_id != "grey")
all_uniprot <- unique(wgcna_mem$protein_id)

suppressMessages({
  sym_map <- clusterProfiler::bitr(
    all_uniprot,
    fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db
  )
})
universe_syms <- unique(sym_map$SYMBOL)
cat(sprintf("ORA universe (full proteome, SYMBOL): %d\n", length(universe_syms)))
uniprot_to_syms <- function(ids) sym_map$SYMBOL[sym_map$UNIPROT %in% ids]

pw <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = FALSE, exclude_variants = TRUE
)

modules <- sort(unique(wgcna_non_grey$group_id))
ora_rows <- list()
for (m in modules) {
  genes <- intersect(uniprot_to_syms(wgcna_non_grey$protein_id[wgcna_non_grey$group_id == m]), universe_syms)
  if (length(genes) < 3) next
  res <- run_ora_deduplicated(
    genes = genes, universe = universe_syms, pathways = pw,
    jaccard_cutoff = 0.5, min_size = 15, max_size = 500, padj_cutoff = 0.05
  )
  if (nrow(res) > 0) ora_rows[[m]] <- res |> mutate(group = m)
}

all_ora <- if (length(ora_rows)) {
  bind_rows(ora_rows) |>
    dplyr::select(group, pathway, database, padj, pval, overlap, size, odds_ratio)
} else {
  tibble(
    group = character(), pathway = character(), database = character(),
    padj = numeric(), pval = numeric(), overlap = integer(),
    size = integer(), odds_ratio = numeric()
  )
}
write.csv(all_ora, file.path(DAT_DIR, "module_ora.csv"), row.names = FALSE)
cat(sprintf("Modules with enrichment: %d / %d\n", length(ora_rows), length(modules)))

top_df <- all_ora |>
  mutate(fill = group) |>
  group_by(group) |>
  slice_min(padj, n = 5, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    log10_padj = -log10(padj),
    clean_name = stringr::str_trunc(clean_pathway_name(pathway), 40, ellipsis = "…"),
    database = if_else(database %in% names(DB_COLORS), database, "GO:BP"),
    yk = reorder_within(clean_name, log10_padj, group)
  )

p_b <- ggplot(top_df, aes(log10_padj, yk)) +
  geom_col(aes(fill = fill, colour = database), linewidth = 0.6, width = 0.82) +
  shadowtext::geom_shadowtext(aes(label = clean_name),
    x = 0.04, hjust = 0, size = 2.1, fontface = "bold",
    colour = "grey10", bg.colour = "white", bg.r = 0.12
  ) +
  facet_wrap(~group, scales = "free", ncol = 3) +
  scale_fill_identity() +
  scale_colour_manual(values = DB_COLORS, name = "Source database") +
  scale_y_reordered() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04)), breaks = scales::breaks_pretty(3)) +
  guides(colour = guide_legend(override.aes = list(fill = "white"))) +
  labs(
    title = "Module biology: pathway over-representation by WGCNA module",
    subtitle = paste(strwrap(sprintf(
      "Hypergeometric ORA with Jaccard de-duplication against the full measured proteome (%d proteins); top 5 pathways per module, adj. p < 0.05. Bar = module colour, outline = source database.",
      length(universe_syms)
    ), width = 125), collapse = "\n"),
    x = expression(-log[10] ~ adj. ~ p), y = NULL, tag = "B"
  ) +
  FIG_THEME +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", size = FIG_AXIS_TEXT),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle = element_text(face = "italic", size = FIG_SUBTITLE_SIZE - 0.5, color = "grey30"),
    plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE)
  )

dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)
save_panel(p_b, file.path(RPT_DIR, "supp", "module_ora_bars"), 250, 230)
cat("F04 module ORA (supplement) done.\n")
