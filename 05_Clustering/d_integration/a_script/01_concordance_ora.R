# Cross-engine concordance (WGCNA vs Mfuzz) and ORA per group.
# Outputs: concordance.csv, concordance_heatmap.pdf, ora_results.csv,
#          ora_dotplots.pdf

pacman::p_load(
  here, dplyr, tidyr, ggplot2,
  clusterProfiler, org.Hs.eg.db, msigdbr
)

source(here("05_Clustering", "functions", "inputs.R"))
i <- load_clustering_inputs()
background_uni <- i$pi_set

w <- read.csv(here(
  "05_Clustering", "a_wgcna_paired", "c_data", "membership.csv"
))
m <- read.csv(here(
  "05_Clustering", "b_mfuzz_gap", "c_data", "membership.csv"
))

# concordance

shared <- inner_join(
  dplyr::select(w, protein_id, wgcna = group_id),
  dplyr::select(m, protein_id, mfuzz = group_id),
  by = "protein_id"
)

ari <- tryCatch(
  {
    ct <- table(shared$wgcna, shared$mfuzz)
    n <- sum(ct)
    sum_row <- rowSums(ct)
    sum_col <- colSums(ct)
    nc2 <- choose(n, 2)
    a <- sum(choose(ct, 2))
    b <- sum(choose(sum_row, 2))
    d <- sum(choose(sum_col, 2))
    (a - b * d / nc2) / ((b + d) / 2 - b * d / nc2)
  },
  error = function(e) NA_real_
)

message(sprintf("Adjusted Rand Index (WGCNA vs Mfuzz): %.4f", ari))

shared_ids <- intersect(w$protein_id, m$protein_id)
w_shared <- w[w$protein_id %in% shared_ids, ]
m_shared <- m[m$protein_id %in% shared_ids, ]

jac_df <- tidyr::crossing(
  wgcna_group = unique(w_shared$group_id),
  mfuzz_group = as.character(unique(m_shared$group_id))
) |>
  mutate(
    n_overlap = mapply(function(wg, mg) {
      wset <- w_shared$protein_id[w_shared$group_id == wg]
      mset <- m_shared$protein_id[
        as.character(m_shared$group_id) == mg
      ]
      length(intersect(wset, mset))
    }, wgcna_group, mfuzz_group),
    n_union = mapply(function(wg, mg) {
      wset <- w_shared$protein_id[w_shared$group_id == wg]
      mset <- m_shared$protein_id[
        as.character(m_shared$group_id) == mg
      ]
      length(union(wset, mset))
    }, wgcna_group, mfuzz_group),
    jaccard = ifelse(n_union > 0, n_overlap / n_union, 0)
  ) |>
  dplyr::select(-n_union)

out_conc <- bind_rows(
  data.frame(
    wgcna_group = "ARI", mfuzz_group = "overall",
    n_overlap = NA_integer_, jaccard = ari
  ),
  jac_df
)
write.csv(
  out_conc,
  here("05_Clustering", "d_integration", "c_data", "concordance.csv"),
  row.names = FALSE
)

heat_p <- ggplot(jac_df, aes(
  x = factor(mfuzz_group), y = wgcna_group, fill = jaccard
)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(jaccard, 2)), size = 3) +
  scale_fill_gradient(low = "white", high = "#2166ac", limits = c(0, 1)) +
  labs(
    x = "Mfuzz cluster", y = "WGCNA module",
    fill = "Jaccard", title = "WGCNA vs Mfuzz overlap"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

heatmap_pdf <- here(
  "05_Clustering", "d_integration", "b_reports", "concordance_heatmap.pdf"
)
pdf(heatmap_pdf, width = 6, height = 4)
print(heat_p)
dev.off()

# ora

reactome_t2g <- msigdbr(
  species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"
) |>
  dplyr::select(gs_name, ncbi_gene) |>
  dplyr::mutate(ncbi_gene = as.character(ncbi_gene)) |>
  dplyr::distinct()

bg_map <- bitr(
  background_uni,
  fromType = "UNIPROT", toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
bg_entrez <- unique(bg_map$ENTREZID)
map_rate <- nrow(bg_map) / length(background_uni)
message(sprintf(
  "Background UniProt->ENTREZ mapping rate: %.1f%% (%d/%d)",
  100 * map_rate, nrow(bg_map), length(background_uni)
))

wgcna_groups <- setdiff(unique(w$group_id), "grey")
mfuzz_groups <- unique(as.character(m$group_id))

run_ora <- function(proteins_uni, grp_entrez_all) {
  grp_map <- bitr(
    proteins_uni,
    fromType = "UNIPROT", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  grp_entrez <- unique(grp_map$ENTREZID)
  if (length(grp_entrez) < 2) {
    return(list(go = list(), reactome = NULL))
  }

  go_res <- lapply(c("BP", "CC", "MF"), function(ont) {
    tryCatch(
      enrichGO(
        gene = grp_entrez,
        universe = grp_entrez_all,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        readable = FALSE
      ),
      error = function(e) NULL
    )
  })
  names(go_res) <- c("GO_BP", "GO_CC", "GO_MF")

  react_res <- tryCatch(
    enricher(
      gene = grp_entrez,
      universe = grp_entrez_all,
      TERM2GENE = reactome_t2g,
      pAdjustMethod = "BH"
    ),
    error = function(e) NULL
  )

  list(go = go_res, reactome = react_res)
}

ora_to_df <- function(res_obj, engine, group_id, term_source) {
  if (is.null(res_obj)) {
    return(NULL)
  }
  df <- as.data.frame(res_obj)
  if (nrow(df) == 0) {
    return(NULL)
  }
  df$engine <- engine
  df$group_id <- as.character(group_id)
  df$source <- term_source
  df[, c(
    "engine", "group_id", "source", "ID", "Description",
    "GeneRatio", "pvalue", "p.adjust", "Count"
  )]
}

n_slots <- (length(wgcna_groups) + length(mfuzz_groups)) * 4
collect <- vector("list", n_slots)
idx <- 0L

for (grp in wgcna_groups) {
  proteins <- w$protein_id[w$group_id == grp]
  res <- run_ora(proteins, bg_entrez)
  for (ont in names(res$go)) {
    idx <- idx + 1L
    collect[[idx]] <- ora_to_df(res$go[[ont]], "wgcna", grp, ont)
  }
  idx <- idx + 1L
  collect[[idx]] <- ora_to_df(res$reactome, "wgcna", grp, "Reactome")
}

for (grp in mfuzz_groups) {
  proteins <- m$protein_id[m$group_id == as.integer(grp)]
  res <- run_ora(proteins, bg_entrez)
  for (ont in names(res$go)) {
    idx <- idx + 1L
    collect[[idx]] <- ora_to_df(res$go[[ont]], "mfuzz", grp, ont)
  }
  idx <- idx + 1L
  collect[[idx]] <- ora_to_df(res$reactome, "mfuzz", grp, "Reactome")
}

ora_df <- bind_rows(Filter(Negate(is.null), collect[seq_len(idx)]))

if (nrow(ora_df) == 0) {
  ora_df <- data.frame(
    engine = character(), group_id = character(), source = character(),
    ID = character(), Description = character(),
    GeneRatio = character(),
    pvalue = numeric(), p.adjust = numeric(), Count = integer()
  )
  message("No ORA results — writing empty ora_results.csv")
}

write.csv(
  ora_df,
  here("05_Clustering", "d_integration", "c_data", "ora_results.csv"),
  row.names = FALSE
)

sig_groups <- if (nrow(ora_df) > 0) {
  unique(paste(
    ora_df$engine[ora_df$p.adjust < 0.05],
    ora_df$group_id[ora_df$p.adjust < 0.05]
  ))
} else {
  character()
}
message(sprintf(
  "Groups with significant enrichment (p.adjust<0.05): %d",
  length(sig_groups)
))

if (nrow(ora_df) > 0) {
  top_terms <- ora_df |>
    filter(p.adjust < 0.05) |>
    group_by(engine, group_id, source) |>
    slice_min(order_by = p.adjust, n = 8, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      neg_log10_padj = -log10(p.adjust),
      label = paste0(engine, "_", group_id, " (", source, ")")
    )

  if (nrow(top_terms) > 0) {
    dot_p <- ggplot(top_terms, aes(
      x = neg_log10_padj,
      y = reorder(Description, neg_log10_padj),
      size = Count, color = source
    )) +
      geom_point() +
      facet_wrap(~label, scales = "free_y", ncol = 2) +
      labs(
        x = "-log10(p.adjust)", y = NULL, size = "Gene count",
        color = "Source", title = "Top ORA terms per group"
      ) +
      theme_bw(base_size = 9) +
      theme(
        strip.text = element_text(size = 7),
        axis.text.y = element_text(size = 6)
      )

    n_facets <- length(unique(top_terms$label))
    ph <- max(4, ceiling(n_facets / 2) * 3)

    pdf(
      here(
        "05_Clustering", "d_integration", "b_reports", "ora_dotplots.pdf"
      ),
      width = 12, height = ph
    )
    print(dot_p)
    dev.off()
  } else {
    pdf(
      here(
        "05_Clustering", "d_integration", "b_reports", "ora_dotplots.pdf"
      ),
      width = 6, height = 4
    )
    plot.new()
    text(0.5, 0.5, "No significant ORA terms (p.adjust < 0.05)")
    dev.off()
  }
} else {
  pdf(
    here(
      "05_Clustering", "d_integration", "b_reports", "ora_dotplots.pdf"
    ),
    width = 6, height = 4
  )
  plot.new()
  text(0.5, 0.5, "No ORA results")
  dev.off()
}

message("Done.")
