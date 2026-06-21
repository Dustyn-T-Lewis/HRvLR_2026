# F07 Panel E: Unified ORA-DEP Chord Diagram -- Training
# Single chord connecting ORA-enriched pathways to "Training" DEPs
# (Training_HR + Training_LR + Training_Interaction, Pi < 0.05 in any).
#
# Ring design (outside -> inside):
#   Track 1: contrast identity color (proteins) | pathway identity color (pathways)
#   Track 2: gene symbol labels     | pathway name (bold white on colored bg)
#   Track 3: logFC bar (proteins)   | padj gradient: UP=red, DOWN=blue (pathways)
#   Ribbons: protein Track 3 -> pathway Track 3, colored by pathway (alpha=0.25)
#
# ORA: Jaccard 0.3 dedup, rrvgo on GO:BP. All surviving pathways shown.
# Direction assigned per pathway via sign(mean logFC of overlap genes).
# Legend column on right side of panel.
# ---------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F07/a_script/style.R")

library(tidyverse)
library(fgsea)
library(circlize)

set.seed(42)

# -- Figure-specific constants ------------------------------------------------

FIGURE_LABEL       <- "Training"
RELEVANT_CONTRASTS <- c("Training_HR", "Training_LR", "Training_Interaction")

RPT <- "04_Figures/F07/b_reports"
DAT <- "04_Figures/F07/c_data"
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

PANEL_W <- 210
PANEL_H <- 180
DPI     <- 300

PW_COLORS_10 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                   "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#FB9A99")

# -- Step 1: Load DEP results & build pool -----------------------------------

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
all_genes <- unique(dep_df$gene)

pi_cols   <- paste0("pi_score_", RELEVANT_CONTRASTS)
logfc_cols <- paste0("logFC_", RELEVANT_CONTRASTS)

# DEP membership per relevant contrast
dep_pool <- dep_df %>%
  mutate(across(all_of(pi_cols), ~ !is.na(.) & . < 0.05,
                .names = "sig_{.col}")) %>%
  filter(if_any(starts_with("sig_pi_score_"), ~ .))

sig_cols <- paste0("sig_pi_score_", RELEVANT_CONTRASTS)
names(dep_pool)[match(sig_cols, names(dep_pool))] <- paste0("sig_", RELEVANT_CONTRASTS)

message(sprintf("DEP pool (%s): %d proteins", FIGURE_LABEL, nrow(dep_pool)))
for (ctr in RELEVANT_CONTRASTS) {
  n <- sum(dep_pool[[paste0("sig_", ctr)]])
  message(sprintf("  %s: %d", ctr, n))
}

# Primary contrast = lowest pi_score among relevant contrasts
primary_df <- dep_pool %>%
  select(gene, all_of(pi_cols)) %>%
  pivot_longer(-gene, names_to = "pi_col", values_to = "pi_val") %>%
  mutate(contrast = gsub("pi_score_", "", pi_col)) %>%
  filter(!is.na(pi_val), contrast %in% RELEVANT_CONTRASTS) %>%
  group_by(gene) %>%
  slice_min(pi_val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, primary_contrast = contrast)

# Display logFC = logFC from primary contrast
display_lfc <- dep_pool %>%
  select(gene, all_of(logfc_cols)) %>%
  left_join(primary_df, by = "gene") %>%
  pivot_longer(all_of(logfc_cols), names_to = "lfc_col", values_to = "lfc_val") %>%
  mutate(lfc_contrast = gsub("logFC_", "", lfc_col)) %>%
  filter(lfc_contrast == primary_contrast) %>%
  select(gene, primary_contrast, display_logFC = lfc_val)

message(sprintf("Primary contrast distribution:"))
print(table(display_lfc$primary_contrast))

# -- Step 2: Membership-based pathway assignment (avoids ORA significance) -----

pw_collection <- build_pathway_collection(min_size = 10, max_size = 500)
pool_genes <- display_lfc$gene

mem_result <- assign_pathways_membership(
  fg_genes     = pool_genes,
  universe     = all_genes,
  pathways     = pw_collection,
  max_pathways = 12, min_overlap = 2,
  jaccard_cutoff = 0.3
)

message(sprintf("Membership assignment: %d / %d mapped (%.0f%%)",
                mem_result$n_mapped, length(pool_genes), mem_result$pct_mapped))

# Filter to non-Other pathways
ora_final <- mem_result$ora %>%
  filter(pathway != "Other", !is.na(pathway))

if (nrow(ora_final) == 0) {
  message("No pathway assignments -- stopping")
  quit(save = "no", status = 0)
}

# Clean pathway names for display
ora_final$pathway_label <- gsub("^HALLMARK_|^REACTOME_|^KEGG_|^GOBP_|^WP_", "",
                                  ora_final$pathway_id %||% ora_final$pathway)
ora_final$pathway_label <- gsub("_", " ", ora_final$pathway_label)
ora_final$pathway_label <- tools::toTitleCase(tolower(ora_final$pathway_label))

# -- Step 3: Build chord data -----------------------------------------------

# Use gene_map to build links
chord_links <- list()
mapped_genes <- mem_result$gene_map %>% filter(pathway != "Other")
pw_names <- unique(mapped_genes$pathway)
for (i in seq_along(pw_names)) {
  pw <- pw_names[i]
  genes <- mapped_genes %>% filter(pathway == pw) %>% pull(gene)
  if (length(genes) == 0) next
  # Look up pathway info from ora results
  pw_info <- mem_result$ora %>% filter(pathway == pw)
  chord_links[[i]] <- tibble(
    gene          = genes,
    pathway       = pw,
    pathway_padj  = if (nrow(pw_info) > 0) pw_info$pval[1] else 1,
    pathway_db    = if (nrow(pw_info) > 0) pw_info$database[1] else "Other",
    pathway_size  = length(genes)
  )
}
link_df <- bind_rows(chord_links)

shorten_pw <- function(x, max_chars = 28) {
  abbrev <- c(
    "Degradation Of The Extracellular Matrix" = "ECM Degradation",
    "Extracellular Matrix Organization"       = "ECM Organization",
    "Collagen Biosynthesis And Modifying Enzymes" = "Collagen Biosynthesis",
    "Integrin Cell Surface Interactions"       = "Integrin Interactions",
    "Epithelial Mesenchymal Transition"        = "EMT",
    "Interleukin 4 And Interleukin 13 Signaling" = "IL-4/IL-13 Signaling",
    "Striated Muscle Contraction"              = "Striated Muscle",
    "Muscle Contraction"                       = "Muscle Contraction",
    "Reference Translation Initiation"         = "Translation Init.",
    "Actin Filament Based Movement"            = "Actin Movement",
    "Muscle Cell Development"                  = "Muscle Cell Dev.",
    "Sarcomere Organization"                   = "Sarcomere Org.",
    "Actomyosin Structure Organization"        = "Actomyosin Org.",
    "Formation Of The Dystrophin Glycoprotein Complex Dgc" = "Dystrophin Complex"
  )
  x <- ifelse(x %in% names(abbrev), abbrev[x], x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 2), ".."), x)
}

link_df <- link_df %>%
  mutate(
    pathway_full  = clean_pathway_name(pathway),
    pathway_label = shorten_pw(clean_pathway_name(pathway))
  )

link_df <- link_df %>%
  left_join(display_lfc, by = "gene")

proteins <- sort(unique(link_df$gene))
n_prot   <- length(proteins)
message(sprintf("Proteins in chord: %d (from %d DEPs)", n_prot, nrow(dep_pool)))

n_missing <- sum(is.na(link_df$primary_contrast))
if (n_missing > 0) {
  warning(sprintf("%d proteins in ORA overlap have no primary contrast -- dropping", n_missing))
  link_df <- link_df %>% filter(!is.na(primary_contrast))
}

# -- Step 4: Direction assignment per pathway --------------------------------

pw_direction <- link_df %>%
  group_by(pathway, pathway_full, pathway_label, pathway_padj, pathway_db,
           pathway_size) %>%
  summarize(
    mean_logFC = mean(display_logFC, na.rm = TRUE),
    n_overlap  = n_distinct(gene),
    .groups    = "drop"
  ) %>%
  mutate(direction = ifelse(mean_logFC > 0, "Up", "Down")) %>%
  arrange(pathway_padj) %>%
  group_by(pathway_full) %>%
  slice_min(pathway_padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pathway_label = make.unique(pathway_label, sep = " ")) %>%
  arrange(pathway_padj)

message("\nPathway direction assignments:")
for (i in seq_len(nrow(pw_direction))) {
  cat(sprintf("  %2d. [%s] %-35s padj=%.1e  mean_lfc=%+.2f  %s  n=%d\n",
    i, pw_direction$pathway_db[i],
    pw_direction$pathway_label[i],
    pw_direction$pathway_padj[i],
    pw_direction$mean_logFC[i],
    pw_direction$direction[i],
    pw_direction$n_overlap[i]))
}

# -- Step 5: Order proteins & pathways ---------------------------------------

pw_direction <- pw_direction %>% arrange(direction, pathway_padj)
pathways_ordered <- pw_direction$pathway_label

gene_primary_pw <- link_df %>%
  group_by(gene) %>%
  slice_min(pathway_padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, primary_pw = pathway_label)

lfc_lookup <- setNames(display_lfc$display_logFC, display_lfc$gene)
ctr_lookup <- setNames(display_lfc$primary_contrast, display_lfc$gene)

proteins_ordered <- gene_primary_pw %>%
  mutate(
    pw_rank  = match(primary_pw, pathways_ordered),
    contrast = ctr_lookup[gene],
    abs_lfc  = abs(lfc_lookup[gene])
  ) %>%
  arrange(pw_rank, contrast, desc(abs_lfc)) %>%
  pull(gene) %>%
  unique()

# -- Step 6: Save chord data CSVs -------------------------------------------

write_csv(
  link_df %>% select(gene, pathway_label, pathway_padj, pathway_db,
                      primary_contrast, display_logFC),
  file.path(DAT, "panel_E", "chord_training_combined.csv")
)

write_csv(
  pw_direction %>% select(pathway_label, pathway_padj, pathway_db,
                            n_overlap, mean_logFC, direction),
  file.path(DAT, "panel_E", "chord_training_pathways.csv")
)
message("\nChord data saved to c_data/panel_E/")

# -- Step 7: Render chord diagram --------------------------------------------

logfc_color <- function(lfc) {
  lfc_clamped <- pmin(pmax(lfc, -2), 2)
  colorRamp2(c(-2, 0, 2), c("#4393C3", "white", "#D6604D"))(lfc_clamped)
}

direction_color <- function(mean_lfc) {
  lfc_clamped <- pmin(pmax(mean_lfc, -2), 2)
  colorRamp2(c(-2, 0, 2), c("#0D47A1", "#F5F5F5", "#B71C1C"))(lfc_clamped)
}

n_pw <- nrow(pw_direction)
pw_col_map <- setNames(PW_COLORS_10[seq_len(n_pw)], pathways_ordered)

pw_dir_col <- setNames(
  vapply(pw_direction$mean_logFC, function(m) direction_color(m), character(1)),
  pathways_ordered
)

pw_sizes <- setNames(
  pmax(pw_direction$n_overlap, 5),
  pathways_ordered
)
all_sectors <- c(proteins_ordered, pathways_ordered)
sector_sizes <- setNames(
  c(rep(3, n_prot), pw_sizes[pathways_ordered]),
  all_sectors
)

prot_gap <- max(0.2, 15 / n_prot)
gap_vec <- c(
  rep(prot_gap, n_prot - 1), 12,
  rep(2, max(n_pw - 1, 0)), 12
)

gene_cex <- max(0.28, min(0.55, 14 / n_prot))
pw_cex   <- 0.50

ctr_col_lookup <- CONTRAST_COLORS[ctr_lookup[proteins_ordered]]

# --- Render function ---
render_chord <- function(filepath, ext) {
  if (ext == "pdf") {
    pdf(filepath, width = PANEL_W / 25.4, height = PANEL_H / 25.4)
  } else {
    png(filepath, width = PANEL_W, height = PANEL_H, units = "mm", res = DPI,
        bg = "white")
  }

  layout(matrix(c(1, 2), ncol = 2), widths = c(0.72, 0.28))

  # --- Panel 1: Chord diagram ---
  par(mar = c(0.5, 0.5, 2.5, 0.5))

  circos.clear()
  circos.par(
    gap.after      = gap_vec,
    start.degree   = 90,
    clock.wise     = TRUE,
    cell.padding   = c(0, 0, 0, 0),
    track.margin   = c(0.005, 0.005)
  )
  circos.initialize(
    factors = factor(all_sectors, levels = all_sectors),
    xlim    = cbind(rep(0, length(all_sectors)), sector_sizes)
  )

  # TRACK 1 (outermost)
  circos.track(ylim = c(0, 1), track.height = 0.14, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index
      if (s %in% proteins_ordered) {
        circos.points(CELL_META$xcenter, 0.12, pch = 15,
                      col = ctr_col_lookup[s], cex = gene_cex * 2.0)
        circos.text(CELL_META$xcenter, 0.55, s,
                    cex = gene_cex, facing = "clockwise", niceFacing = TRUE,
                    font = 3)
      } else if (s %in% pathways_ordered) {
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = pw_dir_col[s], border = "black", lwd = 0.5)
        arc_w <- CELL_META$cell.xlim[2] - CELL_META$cell.xlim[1]
        label_len <- nchar(s)
        dyn_cex <- min(pw_cex, max(0.18, pw_cex * arc_w / max(label_len * 0.7, 8)))
        pw_idx <- match(s, pathways_ordered)
        mean_lfc <- pw_direction$mean_logFC[pw_idx]
        txt_col <- if (abs(mean_lfc) > 0.3) "white" else "grey20"
        circos.text(CELL_META$xcenter, 0.5, s,
                    cex = dyn_cex, font = 2, col = txt_col,
                    facing = "bending.inside", niceFacing = TRUE)
      }
    })

  # TRACK 2: contrast strip (proteins)
  circos.track(ylim = c(0, 1), track.height = 0.040, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index
      if (s %in% proteins_ordered) {
        col <- ctr_col_lookup[s]
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = col, border = NA)
      }
    })

  # TRACK 3 (innermost): logFC bar (proteins) | solid unique color (pathways)
  circos.track(ylim = c(0, 1), track.height = 0.055, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index
      if (s %in% proteins_ordered) {
        col <- logfc_color(lfc_lookup[s])
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = col, border = NA)
      } else if (s %in% pathways_ordered) {
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = pw_col_map[s], border = NA)
      }
    })

  # RIBBONS
  link_unique <- link_df %>%
    filter(gene %in% proteins_ordered, pathway_label %in% pathways_ordered) %>%
    distinct(gene, pathway_label)

  for (i in seq_len(nrow(link_unique))) {
    gene_i <- link_unique$gene[i]
    pw_i   <- link_unique$pathway_label[i]

    pw_genes <- link_unique %>%
      filter(pathway_label == pw_i) %>%
      pull(gene)
    pw_genes <- pw_genes[pw_genes %in% proteins_ordered]
    idx <- which(pw_genes == gene_i)
    if (length(idx) == 0) next

    pw_xlim <- get.cell.meta.data("cell.xlim", sector.index = pw_i,
                                   track.index = 3)
    slot_w <- (pw_xlim[2] - pw_xlim[1]) / length(pw_genes)

    ribbon_col <- adjustcolor(pw_col_map[pw_i], alpha.f = 0.30)

    tryCatch(
      circos.link(
        gene_i, c(0, 1),
        pw_i, c(pw_xlim[1] + (idx - 1) * slot_w,
                pw_xlim[1] + idx * slot_w),
        col    = ribbon_col,
        border = NA,
        rou1   = get.cell.meta.data("cell.bottom.radius",
                                     sector.index = gene_i, track.index = 3),
        rou2   = get.cell.meta.data("cell.bottom.radius",
                                     sector.index = pw_i, track.index = 3)
      ),
      error = function(e) NULL
    )
  }

  title(main = paste(FIGURE_LABEL, "DEP - Pathway Chord"),
        cex.main = 1.3, font.main = 2, line = 0.5)
  circos.clear()

  # --- Panel 2: Legends ---
  par(mar = c(2, 0.5, 2.5, 1))
  plot.new()

  y_cursor <- 0.95

  draw_legend <- function(y_start, title, labels, colors, cex_title = 0.75,
                           cex_labels = 0.65, box_h = 0.022, spacing = 0.028) {
    text(0.05, y_start, title, adj = c(0, 1), font = 2, cex = cex_title)
    y <- y_start - 0.03
    for (i in seq_along(labels)) {
      rect(0.05, y - box_h, 0.15, y, col = colors[i], border = "grey50",
           lwd = 0.3)
      text(0.18, y - box_h / 2, labels[i], adj = c(0, 0.5), cex = cex_labels)
      y <- y - spacing
    }
    return(y - 0.01)
  }

  # 1. Primary Contrast
  y_cursor <- draw_legend(
    y_cursor, "Primary Contrast",
    labels = gsub("_", " ", RELEVANT_CONTRASTS),
    colors = CONTRAST_COLORS[RELEVANT_CONTRASTS]
  )

  # 2. Protein logFC gradient bar
  n_steps <- 50
  x_left <- 0.05; x_right <- 0.55
  bar_h <- 0.018

  text(0.05, y_cursor, "Protein logFC", adj = c(0, 1), font = 2, cex = 0.75)
  y_cursor <- y_cursor - 0.03
  lfc_seq <- seq(-2, 2, length.out = n_steps)
  lfc_cols <- logfc_color(lfc_seq)
  for (j in seq_len(n_steps)) {
    x0 <- x_left + (j - 1) / n_steps * (x_right - x_left)
    x1 <- x_left + j / n_steps * (x_right - x_left)
    rect(x0, y_cursor - bar_h, x1, y_cursor, col = lfc_cols[j], border = NA)
  }
  rect(x_left, y_cursor - bar_h, x_right, y_cursor, col = NA, border = "grey50",
       lwd = 0.3)
  text(x_left, y_cursor - bar_h - 0.012, "-2", cex = 0.55, adj = c(0, 1))
  text((x_left + x_right) / 2, y_cursor - bar_h - 0.012, "0",
       cex = 0.55, adj = c(0.5, 1))
  text(x_right, y_cursor - bar_h - 0.012, "+2", cex = 0.55, adj = c(1, 1))
  y_cursor <- y_cursor - bar_h - 0.04

  # 3. Pathway Direction gradient bar
  text(0.05, y_cursor, "Pathway Direction", adj = c(0, 1), font = 2, cex = 0.75)
  y_cursor <- y_cursor - 0.03
  dir_seq <- seq(-2, 2, length.out = n_steps)
  dir_cols <- direction_color(dir_seq)
  for (j in seq_len(n_steps)) {
    x0 <- x_left + (j - 1) / n_steps * (x_right - x_left)
    x1 <- x_left + j / n_steps * (x_right - x_left)
    rect(x0, y_cursor - bar_h, x1, y_cursor, col = dir_cols[j], border = NA)
  }
  rect(x_left, y_cursor - bar_h, x_right, y_cursor, col = NA,
       border = "grey50", lwd = 0.3)
  text(x_left, y_cursor - bar_h - 0.012, "-2 (Down)", cex = 0.50, adj = c(0, 1))
  text(x_right, y_cursor - bar_h - 0.012, "+2 (Up)", cex = 0.50, adj = c(1, 1))
  y_cursor <- y_cursor - bar_h - 0.04

  # 4. Pathway Identity
  legend_pw_labels <- pw_direction$pathway_label
  y_cursor <- draw_legend(
    y_cursor, "Pathway Identity",
    labels = legend_pw_labels,
    colors = pw_col_map[pathways_ordered],
    cex_labels = 0.50,
    box_h = 0.018,
    spacing = 0.024
  )

  dev.off()
}

# --- Render both formats ---
for (ext in c("pdf", "png")) {
  outfile <- file.path(RPT, paste0("panel_E_ora_chord_MAIN.", ext))
  render_chord(outfile, ext)
  message(sprintf("Rendered: %s", outfile))
}

message(sprintf("\nDone: F07 Panel E -- %s ORA Chord (%d proteins, %d pathways)",
                FIGURE_LABEL, n_prot, n_pw))
