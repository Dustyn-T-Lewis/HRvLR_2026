# Unified pathway enrichment utilities
# MSigDB Hallmark (H), Canonical Pathways (C2:CP), GO:BP (C5:GO:BP); Jaccard dedup per Reimand et al. 2019
#
# Exports:
#   build_pathway_collection()      assemble GO:BP + Reactome + Hallmark + KEGG
#   run_fgsea_deduplicated()        fgsea + collapsePathways for one contrast
#   run_ora_deduplicated()          over-representation analysis with Jaccard dedup
#   classify_database()             database label for plotting

deduplicate_enrichment_flat <- function(results, pathways, jaccard_cutoff = 0.5) {
  if (nrow(results) == 0) {
    return(results)
  }

  results <- results[order(results$padj), ]
  kept_sets <- list()
  keep_mask <- logical(nrow(results))

  for (i in seq_len(nrow(results))) {
    pw_name <- results$pathway[i]
    pw_genes <- pathways[[pw_name]]
    if (is.null(pw_genes)) {
      keep_mask[i] <- TRUE
      next
    }

    is_redundant <- FALSE
    for (j in seq_along(kept_sets)) {
      inter_n <- length(intersect(pw_genes, kept_sets[[j]]))
      union_n <- length(union(pw_genes, kept_sets[[j]]))
      if (union_n > 0 && (inter_n / union_n) > jaccard_cutoff) {
        is_redundant <- TRUE
        break
      }
    }

    if (!is_redundant) {
      keep_mask[i] <- TRUE
      kept_sets[[length(kept_sets) + 1]] <- pw_genes
    }
  }

  results[keep_mask, ]
}

# Database-stratified dedup: within-db first, then merge survivors by padj.
# Falls back to flat dedup if no 'database' column.
deduplicate_enrichment <- function(results, pathways, jaccard_cutoff = 0.5) {
  if (nrow(results) == 0) {
    return(results)
  }

  if (!"database" %in% names(results)) {
    return(deduplicate_enrichment_flat(results, pathways, jaccard_cutoff))
  }

  dbs <- unique(results$database)
  within_dedup <- list()
  for (db in dbs) {
    db_rows <- results[results$database == db, ]
    within_dedup[[db]] <- deduplicate_enrichment_flat(
      db_rows, pathways,
      jaccard_cutoff
    )
  }

  survivors <- do.call(rbind, within_dedup)
  survivors[order(survivors$padj), ]
}


# EnrichmentMap combined similarity (Merico 2010 PMID 21085593; Reimand 2019
# PMID 30664679 Box 1): 0.5 * overlap + 0.5 * jaccard, two sets redundant at
# >= cutoff (0.375 default). This is the display dedup for the F03 rings, kept
# separate from the Jaccard-only cache dedup above, which is frozen because the
# F04/F05 concordance read pivots the cached fgsea_all. Consolidate the two onto
# one coefficient only in a coordinated F03/F04/F05 pass.
em_similarity <- function(a, b) {
  inter <- length(intersect(a, b))
  if (inter == 0L) {
    return(0)
  }
  0.5 * inter / min(length(a), length(b)) +
    0.5 * inter / (length(a) + length(b) - inter)
}

dedup_em_flat <- function(results, pathways, cutoff = 0.375) {
  if (nrow(results) == 0) {
    return(results)
  }
  results <- results[order(results$padj), ]
  kept_sets <- list()
  keep <- logical(nrow(results))
  for (i in seq_len(nrow(results))) {
    genes <- pathways[[results$pathway[i]]]
    if (is.null(genes)) {
      keep[i] <- TRUE
      next
    }
    redundant <- any(vapply(
      kept_sets, function(k) em_similarity(genes, k) >= cutoff, logical(1)
    ))
    if (!redundant) {
      keep[i] <- TRUE
      kept_sets[[length(kept_sets) + 1]] <- genes
    }
  }
  results[keep, ]
}

# Within-database collapse first, then a cross-database pass so the same biology
# in different collections (REACTOME_TCA_CYCLE vs KEGG_CITRATE_CYCLE) reduces to
# one arc.
dedup_em <- function(results, pathways, cutoff = 0.375, cross_db = TRUE) {
  if (nrow(results) == 0) {
    return(results)
  }
  if (!"database" %in% names(results)) {
    return(dedup_em_flat(results, pathways, cutoff))
  }
  within <- lapply(unique(results$database), function(db) {
    dedup_em_flat(results[results$database == db, ], pathways, cutoff)
  })
  survivors <- do.call(rbind, within)
  survivors <- survivors[order(survivors$padj), ]
  if (cross_db && nrow(survivors) > 1) {
    survivors <- dedup_em_flat(survivors, pathways, cutoff)
  }
  survivors
}

# Audit trail: flag each significant pathway as the retained representative or
# redundant, and for redundant ones name the kept pathway it overlaps most and
# the coefficient to it. Reuses dedup_em so the flags cannot drift from the ring.
dedup_report <- function(results, pathways, cutoff = 0.375, cross_db = TRUE) {
  results$dedup_status <- "kept"
  results$merged_into <- NA_character_
  results$overlap_jaccard <- NA_real_
  if (nrow(results) == 0) {
    return(results)
  }
  kept <- dedup_em(results, pathways, cutoff, cross_db)$pathway
  redundant <- setdiff(results$pathway, kept)
  results$dedup_status[results$pathway %in% redundant] <- "redundant"
  for (p in redundant) {
    if (is.null(pathways[[p]])) next
    sims <- vapply(
      kept, function(k) em_similarity(pathways[[p]], pathways[[k]]),
      numeric(1)
    )
    best <- which.max(sims)
    results$merged_into[results$pathway == p] <- kept[best]
    results$overlap_jaccard[results$pathway == p] <- sims[best]
  }
  results[order(results$padj), , drop = FALSE]
}


build_pathway_collection <- function(species = "Homo sapiens",
                                     min_size = 10, max_size = 500,
                                     include_goslim = TRUE,
                                     exclude_variants = FALSE) {
  requireNamespace("msigdbr", quietly = TRUE)

  hallmark <- msigdbr::msigdbr(species = species, collection = "H")
  kegg <- msigdbr::msigdbr(
    species = species, collection = "C2",
    subcollection = "CP:KEGG_MEDICUS"
  )
  reactome <- msigdbr::msigdbr(
    species = species, collection = "C2",
    subcollection = "CP:REACTOME"
  )
  gobp <- msigdbr::msigdbr(
    species = species, collection = "C5",
    subcollection = "GO:BP"
  )
  gocc <- msigdbr::msigdbr(
    species = species, collection = "C5",
    subcollection = "GO:CC"
  )
  gomf <- msigdbr::msigdbr(
    species = species, collection = "C5",
    subcollection = "GO:MF"
  )

  disease_pat <- paste0(
    "DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|",
    "MELANOMA|GLIOMA|HEPATITIS|HIV|INFECTION|VIRAL|",
    "BACTERIAL|PARASIT"
  )
  kegg <- kegg[!grepl(disease_pat, kegg$gs_name, ignore.case = TRUE), ]
  reactome <- reactome[!grepl(disease_pat, reactome$gs_name, ignore.case = TRUE), ]

  if (exclude_variants) {
    kegg <- kegg[!grepl("_VARIANT_", kegg$gs_name), ]
  }

  cols <- c("gs_name", "gene_symbol")
  sets_list <- list(
    hallmark[, cols], kegg[, cols], reactome[, cols],
    gobp[, cols], gocc[, cols], gomf[, cols]
  )
  dbs <- c("H", "KEGG", "Reactome", "GO:BP", "GO:CC", "GO:MF")
  all_sets <- do.call(rbind, sets_list)

  pw_list <- split(all_sets$gene_symbol, all_sets$gs_name)
  pw_list <- lapply(pw_list, unique)

  if (include_goslim) {
    goslim_sets <- build_goslim_gene_sets(
      species = species, min_size = min_size, max_size = max_size
    )
    pw_list <- c(pw_list, goslim_sets)
    dbs <- c(dbs, "GO Slim")
  }

  sizes <- vapply(pw_list, length, integer(1))
  pw_list <- pw_list[sizes >= min_size & sizes <= max_size]

  message(sprintf(
    "Pathway collection: %d sets (%s), size %d-%d",
    length(pw_list), paste(dbs, collapse = " + "),
    min_size, max_size
  ))
  pw_list
}


build_goslim_gene_sets <- function(species = "Homo sapiens",
                                   min_size = 10, max_size = 500) {
  requireNamespace("GO.db", quietly = TRUE)
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  requireNamespace("AnnotationDbi", quietly = TRUE)

  # 62 GO Slim Generic BP terms (from go_slim_categories.R)
  bp_slim <- c(
    "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
    "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
    "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
    "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
    "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
    "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
    "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
    "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
    "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
    "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
    "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
    "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
    "GO:0140014", "GO:1901135"
  )

  # Get all descendant GO terms for each slim term
  offspring <- as.list(GO.db::GOBPOFFSPRING)

  # Map all BP GO terms -> gene symbols via org.Hs.eg.db
  suppressMessages({
    go_genes <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "GO"),
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  })
  go_bp_genes <- go_genes[!is.na(go_genes$ONTOLOGY) & go_genes$ONTOLOGY == "BP", ]
  go_to_symbols <- split(go_bp_genes$SYMBOL, go_bp_genes$GO)

  # Build gene sets: each slim term + all its descendants
  goslim_sets <- list()
  slim_names <- vapply(bp_slim, function(id) {
    tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[id]]),
      error = function(e) NA_character_
    )
  }, character(1))

  for (i in seq_along(bp_slim)) {
    go_id <- bp_slim[i]
    go_term <- slim_names[i]
    if (is.na(go_term)) next

    # Collect genes from this term + all offspring
    all_terms <- go_id
    desc <- offspring[[go_id]]
    if (!is.null(desc)) all_terms <- c(all_terms, desc)

    genes <- unique(unlist(go_to_symbols[intersect(all_terms, names(go_to_symbols))],
      use.names = FALSE
    ))
    genes <- genes[!is.na(genes)]

    if (length(genes) >= min_size && length(genes) <= max_size) {
      set_name <- paste0("GOSLIM_", toupper(gsub(" ", "_", go_term)))
      goslim_sets[[set_name]] <- genes
    }
  }

  message(sprintf(
    "GO Slim: %d/%d terms passed size filter (%d-%d)",
    length(goslim_sets), length(bp_slim), min_size, max_size
  ))
  goslim_sets
}


run_fgsea_deduplicated <- function(ranks, pathways, jaccard_cutoff = 0.5,
                                   nperm = 10000, min_size = 15,
                                   max_size = 500) {
  requireNamespace("fgsea", quietly = TRUE)

  res <- fgsea::fgseaMultilevel(
    pathways    = pathways,
    stats       = ranks,
    minSize     = min_size,
    maxSize     = max_size,
    nPermSimple = nperm,
    eps         = 0
  )
  res <- as.data.frame(res)

  res$database <- classify_database(res$pathway)

  res <- tibble::as_tibble(res)

  keep_cols <- c(
    "pathway", "padj", "NES", "size", "leadingEdge",
    "database", "pval", "ES", "log2err"
  )
  res <- res[, intersect(keep_cols, names(res))]

  sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
  nonsig <- res[is.na(res$padj) | res$padj >= 0.05, ]

  sig_dedup <- deduplicate_enrichment(sig, pathways, jaccard_cutoff)

  n_removed <- nrow(sig) - nrow(sig_dedup)
  pct <- if (nrow(sig) > 0) round(100 * n_removed / nrow(sig), 1) else 0
  message(sprintf(
    "fGSEA dedup: %d sig -> %d kept (removed %d, %.1f%%)",
    nrow(sig), nrow(sig_dedup), n_removed, pct
  ))

  rbind(sig_dedup, nonsig)
}


run_ora_deduplicated <- function(genes, universe, pathways,
                                 jaccard_cutoff = 0.5,
                                 min_size = 10, max_size = 500,
                                 padj_cutoff = 0.05) {
  requireNamespace("fgsea", quietly = TRUE)

  genes <- intersect(genes, universe)

  pw_by_db <- split(names(pathways), classify_database(names(pathways)))
  db_results <- list()
  for (db in names(pw_by_db)) {
    db_pw <- pathways[pw_by_db[[db]]]
    if (length(db_pw) < 2) next
    db_res <- fgsea::fora(
      pathways = db_pw,
      genes    = genes,
      universe = universe,
      minSize  = min_size,
      maxSize  = max_size
    )
    db_res <- as.data.frame(db_res)
    if (nrow(db_res) == 0) next
    db_res$database <- db
    db_results[[db]] <- db_res
  }
  res <- do.call(rbind, db_results)

  N <- length(universe)
  K <- length(genes)
  res$odds_ratio <- vapply(seq_len(nrow(res)), function(i) {
    a <- res$overlap[i]
    b <- K - a
    c <- res$size[i] - a
    d <- N - K - c
    if (b == 0 || c == 0) Inf else (a * d) / (b * c)
  }, numeric(1))

  res <- tibble::as_tibble(res)

  sig <- res[!is.na(res$padj) & res$padj < padj_cutoff, ]
  sig_dedup <- deduplicate_enrichment(sig, pathways, jaccard_cutoff)

  n_removed <- nrow(sig) - nrow(sig_dedup)
  pct <- if (nrow(sig) > 0) round(100 * n_removed / nrow(sig), 1) else 0
  message(sprintf(
    "ORA dedup: %d sig -> %d kept (removed %d, %.1f%%)",
    nrow(sig), nrow(sig_dedup), n_removed, pct
  ))

  sig_dedup
}


classify_database <- function(pathway_names) {
  dplyr::case_when(
    grepl("^HALLMARK_", pathway_names) ~ "Hallmark",
    grepl("^REACTOME_", pathway_names) ~ "Reactome",
    grepl("^KEGG_MEDICUS_", pathway_names) ~ "KEGG",
    grepl("^KEGG_", pathway_names) ~ "KEGG",
    grepl("^GOSLIM_", pathway_names) ~ "GO Slim",
    grepl("^GOBP_", pathway_names) ~ "GO:BP",
    grepl("^GOCC_", pathway_names) ~ "GO:CC",
    grepl("^GOMF_", pathway_names) ~ "GO:MF",
    TRUE ~ "Other"
  )
}
