# F03 Panel D: ORA by module and cluster ----
setwd(rprojroot::find_rstudio_root_file())
if (!exists("wgcna_mem")) source("04_Figures/F03/a_script/HRvLR_F03_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# UniProt -> SYMBOL mapping ----
# Build the full set of UniProt IDs to map: all group proteins + pi_set background
wgcna_non_grey <- wgcna_mem |> filter(group_id != "grey")
mfuzz_labelled <- mfuzz_mem |> mutate(cluster = paste0("c", group_id))

all_uniprot <- unique(c(wgcna_non_grey$protein_id, mfuzz_labelled$protein_id, pi_set))

suppressMessages({
  sym_map <- clusterProfiler::bitr(
    all_uniprot,
    fromType = "UNIPROT",
    toType   = "SYMBOL",
    OrgDb    = org.Hs.eg.db
  )
})
map_rate <- round(100 * nrow(sym_map) / length(all_uniprot), 1)
cat(sprintf("UniProt->SYMBOL mapping: %d/%d (%.1f%%)\n", nrow(sym_map), length(all_uniprot), map_rate))

# Universe: pi_set mapped to SYMBOL (methodologically correct background)
universe_syms <- sym_map$SYMBOL[sym_map$UNIPROT %in% pi_set]
cat(sprintf("Universe size (pi-gated, SYMBOL): %d\n", length(universe_syms)))

# Build pathway collection (gene symbols, no GO Slim)
pw <- build_pathway_collection(
  min_size = 15,
  max_size = 500,
  include_goslim = FALSE,
  exclude_variants = TRUE
)

# Helper: UniProt IDs -> mapped SYMBOL vector
uniprot_to_syms <- function(uniprot_ids) {
  sym_map$SYMBOL[sym_map$UNIPROT %in% uniprot_ids]
}

# Define groups: 5 WGCNA modules + 4 Mfuzz clusters
wgcna_modules <- sort(unique(wgcna_non_grey$group_id))
mfuzz_clusters <- paste0("c", sort(unique(mfuzz_labelled$cluster |>
  stringr::str_remove("^c") |> as.integer())))

groups_list <- c(
  setNames(
    lapply(wgcna_modules, function(m) wgcna_non_grey$protein_id[wgcna_non_grey$group_id == m]),
    paste("wgcna", wgcna_modules)
  ),
  setNames(
    lapply(mfuzz_clusters, function(cl) {
      mfuzz_labelled$protein_id[mfuzz_labelled$cluster == cl]
    }),
    paste("mfuzz", mfuzz_clusters)
  )
)

# Run ORA per group ----
ora_rows <- list()
for (grp in names(groups_list)) {
  genes_sym <- uniprot_to_syms(groups_list[[grp]])
  genes_in_universe <- intersect(genes_sym, universe_syms)

  if (length(genes_in_universe) < 3) {
    message(sprintf("  %s: only %d genes in universe, skipping", grp, length(genes_in_universe)))
    next
  }

  message(sprintf("  Running ORA for %s (%d genes in universe)", grp, length(genes_in_universe)))
  res <- run_ora_deduplicated(
    genes = genes_in_universe,
    universe = universe_syms,
    pathways = pw,
    jaccard_cutoff = 0.5,
    min_size = 15,
    max_size = 500,
    padj_cutoff = 0.05
  )

  if (nrow(res) > 0) {
    engine <- strsplit(grp, " ")[[1]][1]
    group <- strsplit(grp, " ")[[1]][2]
    ora_rows[[grp]] <- res |>
      mutate(engine = engine, group = group, group_label = grp)
  }
}

n_sig_groups <- length(ora_rows)
cat(sprintf("Groups with significant enrichment: %d / %d\n", n_sig_groups, length(groups_list)))

# Combine all results for audit export ----
if (n_sig_groups > 0) {
  all_ora <- bind_rows(ora_rows) |>
    dplyr::select(engine, group, group_label, pathway, database, padj, pval, overlap, size, odds_ratio)
} else {
  all_ora <- tibble(
    engine = character(), group = character(), group_label = character(),
    pathway = character(), database = character(),
    padj = numeric(), pval = numeric(),
    overlap = integer(), size = integer(), odds_ratio = numeric()
  )
}

write.csv(all_ora, file.path(DAT_DIR, "audit_panel_D_ora.csv"), row.names = FALSE)
cat(sprintf("Audit written: %d significant ORA rows across %d groups\n", nrow(all_ora), n_sig_groups))

# Plot ----
if (n_sig_groups > 0) {
  top_df <- all_ora |>
    group_by(group_label) |>
    slice_min(order_by = padj, n = 5, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      log10_padj  = -log10(padj),
      clean_name  = clean_pathway_name(pathway),
      group_label = factor(group_label, levels = names(groups_list))
    )

  # Ensure database column maps to DB_COLORS keys
  db_levels <- names(DB_COLORS)
  top_df <- top_df |>
    mutate(database = if_else(database %in% db_levels, database, "GO:BP"))

  n_facets <- n_sig_groups
  plot_height <- max(150, n_facets * 40)

  p <- ggplot(
    top_df,
    aes(
      x = log10_padj,
      y = reorder_within(clean_name, log10_padj, group_label)
    )
  ) +
    geom_col(
      aes(fill = database),
      colour = "black",
      linewidth = 0.3,
      width = 0.75
    ) +
    scale_y_reordered() +
    scale_fill_manual(values = DB_COLORS, name = "Database") +
    facet_wrap(~group_label, scales = "free_y", ncol = 3) +
    labs(
      x = "-log10 adj. p",
      y = NULL,
      title = "Pathway over-representation by module and cluster",
      subtitle = paste0(
        "Hypergeometric + Jaccard-dedup vs the pi-gated set (298 proteins); ",
        "top 5 per group, adj. p<0.05"
      ),
      tag = "D"
    ) +
    FIG_THEME +
    theme(axis.text.y = element_text(size = 7))
} else {
  # No enrichment — honest empty panel
  p <- ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.5,
      label = "No significant pathways (adj. p < 0.05)\nagainst pi-gated universe (n = 298)",
      size = 4, hjust = 0.5, vjust = 0.5, color = "grey40"
    ) +
    labs(
      title = "Pathway over-representation by module and cluster",
      subtitle = paste0(
        "Hypergeometric + Jaccard-dedup vs the pi-gated set (298 proteins); ",
        "top 5 per group, adj. p<0.05"
      ),
      tag = "D"
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE, color = "grey30"),
      plot.tag      = element_text(face = "bold", size = FIG_TAG_SIZE)
    )

  plot_height <- 80
}

save_panel(p, file.path(RPT_DIR, "panel_d_ora"), 250, plot_height)
cat("F03 Panel D done.\n")
