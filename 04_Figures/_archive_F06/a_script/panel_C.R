# F06 Panel C: Hub Protein Networks (2x2 grid for top 4 modules)
# Nodes = hub genes (kME >= Q90), sized by kME
# Edges = TOM co-expression

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

library(readr)
library(dplyr)
library(tibble)
library(igraph)
library(ggraph)
library(patchwork)

RPT <- "04_Figures/F06/b_reports"
DAT <- "04_Figures/F06/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load data ---
kME           <- readRDS(file.path(DAT, "kME_all.rds"))
module_colors <- readRDS(file.path(DAT, "module_colors.rds"))
module_df     <- read_csv(file.path(DAT, "wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_labels    <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
key_mods      <- readLines(file.path(DAT, "key_modules.txt"))
ann           <- read_csv(file.path(DAT, "imp_annotations.csv"), show_col_types = FALSE)

# Load TOM or adjacency from WGCNA network
net_path <- "05_WGCNA/c_data/wgcna/wgcna_network.rds"
sft_path <- file.path(DAT, "wgcna_sft_summary.csv")

# --- Hub network builder ---
make_hub_network <- function(mod_color) {
  mod_label <- mod_labels$display_label[mod_labels$module_color == mod_color]
  if (length(mod_label) == 0) mod_label <- mod_color

  # Get module genes
  mod_genes <- names(module_colors)[module_colors == mod_color]
  kme_col   <- paste0("kME", mod_color)

  if (!(kme_col %in% colnames(kME))) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = paste(mod_color, "- no kME")))
  }

  # Hub selection: kME >= Q90 within module
  kme_vals <- kME[mod_genes, kme_col, drop = TRUE]
  kme_vals <- kme_vals[!is.na(kme_vals)]
  q90      <- quantile(kme_vals, 0.90)
  hubs     <- names(kme_vals[kme_vals >= q90])

  if (length(hubs) < 3) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = paste(mod_color, "- <3 hubs")))
  }

  # Build adjacency submatrix from correlation
  datExpr <- readRDS(file.path(DAT, "datExpr.rds"))
  hub_expr <- datExpr[, hubs, drop = FALSE]
  cor_mat  <- cor(hub_expr, use = "pairwise.complete.obs")
  diag(cor_mat) <- 0
  cor_mat[cor_mat < 0.3] <- 0  # threshold weak edges

  # Build graph
  g <- graph_from_adjacency_matrix(abs(cor_mat), mode = "undirected",
                                    weighted = TRUE, diag = FALSE)
  g <- delete_edges(g, which(E(g)$weight == 0))

  if (ecount(g) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = paste(mod_color, "- no edges")))
  }

  # Node attributes
  V(g)$kME <- kme_vals[V(g)$name]
  V(g)$gene <- ann$gene[match(V(g)$name, ann$uniprot_id)]
  V(g)$gene[is.na(V(g)$gene)] <- V(g)$name[is.na(V(g)$gene)]

  set.seed(42)
  p <- ggraph(g, layout = "stress") +
    geom_edge_link(aes(alpha = weight), color = "grey60", width = 0.3,
                   show.legend = FALSE) +
    geom_node_point(aes(size = kME), fill = mod_color,
                    shape = 21, color = "black", stroke = 0.3) +
    geom_node_text(aes(label = gene), size = 2, repel = TRUE,
                   max.overlaps = 15) +
    scale_size_continuous(range = c(2, 6), name = "kME") +
    labs(title = paste0(mod_label, " (", mod_color, ", n=", length(hubs), " hubs)")) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
          legend.position = "none")

  p
}

# --- Build 2x2 grid ---
plots <- lapply(key_mods[1:min(4, length(key_mods))], make_hub_network)

pC <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Hub Protein Networks",
    subtitle = "Top hubs (kME \u2265 Q90) | Edge = co-expression (r \u2265 0.3)",
    tag_levels = list(c("C", "", "", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      plot.tag      = element_text(face = "bold", size = 15)
    )
  )

PC_W <- 300; PC_H <- 300
ggsave(file.path(RPT, "panel_C_hub_network_MAIN.pdf"), pC,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_C_hub_network_MAIN.png"), pC,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)
message("F06 Panel C done")
