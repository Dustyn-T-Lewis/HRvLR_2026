################################################################################
#   F09 Panel D: Hub Protein Networks (geom_mark_hull overlay)
#   Key modules, individual output per module. Pathway hulls drawn directly
#   onto the network via ggforce::geom_mark_hull.
#   Hub selection: kME >= module Q90 (data-driven)
#   Pathway DB: run_ora_deduplicated() with full multi-DB collection
#   Layout: stress (graphlayouts, deterministic)
#   Generates: panel_D_hub_{mod}_MAIN.pdf/png, c_data/04_panel_D_*.csv
################################################################################

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F09/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(ggrepel)
  library(WGCNA); library(igraph); library(ggraph)
  library(ggforce); library(concaveman); library(graphlayouts)
  library(tidygraph); library(ggnewscale); library(fgsea); library(colorspace)
})

allowWGCNAThreads()
set.seed(42)

RPT <- "04_Figures/F09/b_reports"
DAT <- "04_Figures/F09/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load data ---
meta       <- read_csv(file.path(PANEL_DATA, "meta.csv"), show_col_types = FALSE)
MEs        <- readRDS(file.path(PANEL_DATA, "MEs.rds"))
kME_all    <- readRDS(file.path(PANEL_DATA, "kME_all.rds"))
datExpr    <- readRDS(file.path(PANEL_DATA, "datExpr.rds"))
module_df  <- read_csv(file.path(WGCNA_DATA, "wgcna_module_assignments.csv"),
                        show_col_types = FALSE)
sft_csv    <- read_csv(file.path(WGCNA_DATA, "wgcna_sft_summary.csv"),
                        show_col_types = FALSE)
NET_POWER  <- sft_csv$selected_power[1]

mod_bio_labels_df  <- read_csv(file.path(PANEL_DATA, "mod_bio_labels.csv"),
                                show_col_types = FALSE)
mod_bio_labels_vec <- setNames(mod_bio_labels_df$bio_label, mod_bio_labels_df$module_color)
display_label_vec  <- setNames(mod_bio_labels_df$display_label, mod_bio_labels_df$module_color)

# Numeric phenotype for GS: group membership (HR=1, LR=0)
meta$group_num <- ifelse(meta$group == "HR", 1, 0)

bg_genes    <- unique(module_df$gene)
KEY_MODULES <- readLines(file.path(PANEL_DATA, "key_modules.txt"))
KEY_MODULES <- KEY_MODULES[nzchar(trimws(KEY_MODULES))]
uid2gene    <- setNames(module_df$gene, module_df$uniprot_id)

# --- Contrasting hull palette (Dark2-derived) ---
HULL_PALETTE <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                  "#66A61E", "#E6AB02", "#A6761D", "#666666")

message("Panel D: hub protein networks...")

# --- Pathway collection (full multi-DB) ---
pw_full <- build_pathway_collection(min_size = 15, max_size = 500)

# --- Hub selection: Q90 per module ---
select_hubs_q90 <- function(mod) {
  mod_prots <- module_df$uniprot_id[module_df$module_color == mod]
  kme_col   <- paste0("kME", mod)
  matched   <- intersect(mod_prots, rownames(kME_all))
  mod_kme   <- setNames(kME_all[matched, kme_col], matched)
  mod_kme   <- mod_kme[!is.na(mod_kme)]
  q90       <- quantile(mod_kme, 0.90)
  names(mod_kme[mod_kme >= q90])
}

# --- Functional group assignment via ORA ---
assign_groups_ora <- function(gene_names, max_groups = 4, min_group_n = 3) {
  clean_pw_name <- function(name) {
    name %>%
      gsub("^HALLMARK_|^GOSLIM_|^GOBP_|^REACTOME_|^KEGG_MEDICUS_", "", .) %>%
      gsub("_", " ", .) %>% str_to_title() %>% str_trunc(35)
  }

  ora_res <- tryCatch(
    run_ora_deduplicated(
      genes    = gene_names,
      universe = bg_genes,
      pathways = pw_full,
      jaccard_cutoff = 0.5,
      min_size = 10, max_size = 500,
      padj_cutoff = 0.05
    ),
    error = function(e) { message("  ORA error: ", e$message); NULL }
  )

  if (is.null(ora_res) || nrow(ora_res) == 0) {
    message("  No significant ORA results")
    return(setNames(rep("Other", length(gene_names)), gene_names))
  }

  ora_res <- ora_res[order(ora_res$padj), ]
  gene_map <- data.frame(gene = character(), pathway = character(),
                         stringsAsFactors = FALSE)
  for (i in seq_len(nrow(ora_res))) {
    hits <- intersect(ora_res$overlapGenes[[i]], gene_names)
    if (length(hits))
      gene_map <- rbind(gene_map, data.frame(gene = hits, pathway = ora_res$pathway[i],
                                             stringsAsFactors = FALSE))
  }
  gene_map <- gene_map[!duplicated(gene_map$gene), ]

  term_counts <- table(gene_map$pathway)
  keep <- names(term_counts[term_counts >= min_group_n])
  keep <- head(keep[order(term_counts[keep], decreasing = TRUE)], max_groups)

  if (length(keep) < 2) {
    message("  Fewer than 2 qualifying groups")
    return(setNames(rep("Other", length(gene_names)), gene_names))
  }

  assignments <- setNames(rep("Other", length(gene_names)), gene_names)
  for (g in gene_names) {
    row <- gene_map[gene_map$gene == g, ]
    if (nrow(row) > 0 && row$pathway[1] %in% keep)
      assignments[g] <- clean_pw_name(row$pathway[1])
  }
  assignments
}

# --- Signed GS (correlation with group membership: HR=1, LR=0) ---
compute_gs_signed <- function() {
  pheno_vec <- meta$group_num
  names(pheno_vec) <- meta$sample_id
  valid_samps <- intersect(meta$sample_id, rownames(datExpr))
  gs <- cor(datExpr[valid_samps, , drop = FALSE], pheno_vec[valid_samps],
            use = "pairwise.complete.obs")
  setNames(gs[, 1], rownames(gs))
}

gs_signed <- compute_gs_signed()

# --- BUILD ONE MODULE ---
build_network_hull <- function(mod) {

  message(sprintf("\n=== %s ===", toupper(mod)))

  hub_ids <- select_hubs_q90(mod)
  hub_genes <- uid2gene[hub_ids]; hub_genes <- hub_genes[!is.na(hub_genes)]
  hub_ids <- hub_ids[hub_ids %in% names(hub_genes)]
  n_mod <- sum(module_df$module_color == mod)
  n_hubs <- length(hub_ids)
  message(sprintf("  %d hubs (Q90 of %d)", n_hubs, n_mod))

  # TOM subnetwork (module-specific)
  cor_fn <- WGCNA::cor
  mod_prots <- intersect(module_df$uniprot_id[module_df$module_color == mod],
                         colnames(datExpr))
  adj_mod <- adjacency(datExpr[, mod_prots], power = NET_POWER, type = "signed hybrid",
                       corFnc = "cor")
  tom_mod <- TOMsimilarity(adj_mod, TOMType = "signed")
  colnames(tom_mod) <- rownames(tom_mod) <- mod_prots
  cor <- stats::cor  # restore

  tom_sub <- tom_mod[hub_ids, hub_ids]
  tom_q90 <- quantile(tom_sub[upper.tri(tom_sub)], 0.90)

  g <- graph_from_adjacency_matrix(tom_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
  E(g)$weight_orig <- E(g)$weight
  g <- delete_edges(g, which(E(g)$weight < tom_q90))
  iso <- which(degree(g) == 0)
  if (length(iso)) g <- delete_vertices(g, iso)

  node_uids  <- V(g)$name
  node_genes <- uid2gene[node_uids]
  node_kme   <- setNames(kME_all[node_uids, paste0("kME", mod)], node_uids)
  node_gs    <- gs_signed[node_uids]; node_gs[is.na(node_gs)] <- 0

  # Labels: top 6 kME + top 3 |GS|
  top_kme <- names(sort(node_kme, decreasing = TRUE))[1:min(6, length(node_kme))]
  top_gs  <- names(sort(abs(node_gs), decreasing = TRUE))[1:min(3, length(node_gs))]
  label_ids <- unique(c(top_kme, top_gs))

  # Functional groups via ORA
  groups <- assign_groups_ora(node_genes)
  node_groups <- groups[node_genes]

  V(g)$gene <- node_genes; V(g)$kME <- node_kme
  V(g)$gs <- node_gs; V(g)$func_grp <- node_groups

  # Stress layout
  set.seed(42)
  lay <- layout_with_stress(g)
  lay_df <- data.frame(x = lay[, 1], y = lay[, 2], name = V(g)$name)

  nd <- data.frame(x = lay_df$x, y = lay_df$y, name = node_uids,
                   gene = node_genes, kME = node_kme, GS = node_gs,
                   func_grp = node_groups, stringsAsFactors = FALSE)

  grp_counts <- table(nd$func_grp)
  nd$n_in_grp <- as.integer(grp_counts[nd$func_grp])

  grp_names <- setdiff(unique(nd$func_grp[nd$func_grp != "Other" & nd$n_in_grp >= 3]), NA)
  hull_colors <- setNames(HULL_PALETTE[seq_along(grp_names)], grp_names)

  tg <- as_tbl_graph(g)

  disp_label <- if (mod %in% names(display_label_vec)) display_label_vec[[mod]]
                else str_to_title(mod)

  hull_nd <- nd %>% filter(func_grp != "Other", n_in_grp >= 3)

  p <- ggraph(tg, layout = "manual", x = lay_df$x, y = lay_df$y)

  # Layer 1: pathway hulls
  if (nrow(hull_nd) > 0 && length(grp_names) > 0) {
    p <- p +
      geom_mark_hull(
        data = hull_nd,
        aes(x = x, y = y, group = func_grp, fill = func_grp),
        concavity = 2, expand = unit(2, "mm"), radius = unit(2, "mm"),
        alpha = 0.15, linewidth = 0.6, show.legend = TRUE,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = hull_colors, name = "Pathway")
  }

  # Layer 2: TOM edges
  p <- p +
    geom_edge_link(aes(width = weight_orig), alpha = 0.55,
                   color = "grey30", show.legend = FALSE) +
    scale_edge_width_continuous(range = c(0.5, 2.0), guide = "none")

  # Layer 3: nodes
  p <- p +
    new_scale_fill() +
    geom_node_point(aes(size = kME, fill = gs), shape = 21,
                    color = "black", stroke = 0.5) +
    scale_size_continuous(range = c(2.0, 7.0), guide = "none") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "GS (Group)",
                         guide = guide_colorbar(barwidth = unit(3, "mm"),
                                                barheight = unit(20, "mm")))

  # Layer 4: gene labels
  p <- p +
    geom_label_repel(
      data = nd[nd$name %in% label_ids, ],
      aes(x = x, y = y, label = gene),
      size = 3.2, fontface = "bold.italic",
      fill = alpha("white", 0.88), color = "grey10",
      label.size = 0.15, label.padding = unit(1.0, "mm"),
      segment.size = 0.25, segment.color = "grey40",
      box.padding = 0.45, point.padding = 0.2,
      max.overlaps = 20, seed = 42, inherit.aes = FALSE
    )

  p <- p +
    labs(title    = disp_label,
         subtitle = sprintf("%d hubs | GS: HR vs LR group membership", n_hubs)) +
    theme_void() +
    theme(
      plot.title      = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle   = element_text(size = 7, hjust = 0.5, color = "grey40"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(4, 4, 4, 4),
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 8),
      legend.text     = element_text(size = 7)
    )

  list(plot = p, node_data = nd)
}

# --- BUILD ALL KEY MODULES ---
results <- lapply(KEY_MODULES, function(mod) build_network_hull(mod))
names(results) <- KEY_MODULES

# --- Export individual panels ---
PD_W <- 170; PD_H <- 170
for (mod in KEY_MODULES) {
  res <- results[[mod]]
  fname <- sprintf("panel_D_hub_%s_MAIN", mod)
  ggsave(file.path(RPT, paste0(fname, ".png")), res$plot,
         width = PD_W, height = PD_H, units = "mm", dpi = 300, limitsize = FALSE)
  ggsave(file.path(RPT, paste0(fname, ".pdf")), res$plot,
         width = PD_W, height = PD_H, units = "mm",
         device = pdf_device, limitsize = FALSE)
  message(sprintf("  Saved %s", fname))
}

# --- Export hub node data ---
all_node_df <- bind_rows(lapply(KEY_MODULES, function(mod) {
  nd <- results[[mod]]$node_data
  if (is.null(nd) || nrow(nd) == 0) return(NULL)
  tibble(uniprot_id = nd$name, module = mod, gene = nd$gene,
         kME = nd$kME, GS = nd$GS, functional_group = nd$func_grp)
}))
write_csv(all_node_df, file.path(DAT, "04_panel_D_hub_network.csv"))

# --- Hub CI computation (stat audit) ---
hub_ci_list <- list()
for (mod in KEY_MODULES) {
  nd <- results[[mod]]$node_data
  if (is.null(nd) || nrow(nd) == 0) next
  me_col <- paste0("ME", mod)
  pheno_vec <- meta$group_num
  names(pheno_vec) <- meta$sample_id
  pheno_vec <- pheno_vec[intersect(names(pheno_vec), rownames(datExpr))]

  for (j in seq_len(nrow(nd))) {
    uid <- nd$name[j]
    if (!(uid %in% colnames(datExpr))) next
    prot_expr <- datExpr[, uid]
    kme_ct <- cor.test(prot_expr, MEs[, me_col], method = "pearson")
    valid <- intersect(names(pheno_vec), names(prot_expr))
    gs_ct <- if (length(valid) >= 4) cor.test(prot_expr[valid], pheno_vec[valid],
                                               method = "pearson") else NULL
    hub_ci_list <- c(hub_ci_list, list(tibble(
      module = mod, uniprot_id = uid, gene = nd$gene[j],
      kME = round(kme_ct$estimate, 4),
      kME_ci_lo = round(kme_ct$conf.int[1], 4),
      kME_ci_hi = round(kme_ct$conf.int[2], 4),
      kME_p = kme_ct$p.value,
      GS = if (!is.null(gs_ct)) round(gs_ct$estimate, 4) else NA_real_,
      GS_ci_lo = if (!is.null(gs_ct)) round(gs_ct$conf.int[1], 4) else NA_real_,
      GS_ci_hi = if (!is.null(gs_ct)) round(gs_ct$conf.int[2], 4) else NA_real_,
      GS_p = if (!is.null(gs_ct)) gs_ct$p.value else NA_real_,
      GS_n = length(valid)
    )))
  }
}
write_csv(bind_rows(hub_ci_list), file.path(DAT, "04_panel_D_hub_CIs.csv"))

message("  Panel D saved (individual modules)")
