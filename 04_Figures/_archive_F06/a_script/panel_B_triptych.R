# F06 Panel B: Per-Module Triptych (z-score heatmap + eigengene + ORA)
# Top 4 modules by LMM significance
# Groups: HR_T1, HR_T2, LR_T1, LR_T2

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)

RPT <- "04_Figures/F06/b_reports"
DAT <- "04_Figures/F06/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load panel-ready data ---
MEs           <- readRDS(file.path(DAT, "MEs.rds"))
datExpr       <- readRDS(file.path(DAT, "datExpr.rds"))
module_colors <- readRDS(file.path(DAT, "module_colors.rds"))
meta          <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)
mod_labels    <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
enrich_df     <- read_csv(file.path(DAT, "wgcna_module_enrichment.csv"), show_col_types = FALSE)
key_mods      <- readLines(file.path(DAT, "key_modules.txt"))

GROUP_ORDER <- c("HR_T1", "HR_T2", "LR_T1", "LR_T2")
GROUP_COLS  <- c(HR_T1 = "#2166AC", HR_T2 = "#67A9CF",
                 LR_T1 = "#B2182B", LR_T2 = "#EF8A62")

# --- Build one triptych row ---
make_triptych_row <- function(mod_color) {
  mod_genes <- names(module_colors)[module_colors == mod_color]
  mod_label <- mod_labels$display_label[mod_labels$module_color == mod_color]
  if (length(mod_label) == 0) mod_label <- mod_color

  # Left: z-score heatmap
  expr_sub <- datExpr[, mod_genes, drop = FALSE]
  z_mat <- scale(expr_sub)
  z_long <- as.data.frame(z_mat) |>
    mutate(sample = rownames(z_mat)) |>
    pivot_longer(-sample, names_to = "gene", values_to = "z") |>
    left_join(meta |> select(Col_ID, Group_Time), by = c("sample" = "Col_ID")) |>
    filter(Group_Time %in% GROUP_ORDER) |>
    mutate(Group_Time = factor(Group_Time, levels = GROUP_ORDER))

  z_summary <- z_long |>
    group_by(Group_Time, gene) |>
    summarise(z = mean(z, na.rm = TRUE), .groups = "drop")

  p_heat <- ggplot(z_summary, aes(x = Group_Time, y = gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "z") +
    labs(title = paste0(mod_label, " (", mod_color, ")"), x = NULL, y = NULL) +
    FIG_THEME +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "none", plot.title = element_text(size = 9))

  # Middle: eigengene trajectory
  me_col <- paste0("ME", mod_color)
  if (me_col %in% colnames(MEs)) {
    me_df <- data.frame(me = MEs[[me_col]], sample = rownames(MEs)) |>
      left_join(meta |> select(Col_ID, Group_Time, Subject_ID),
                by = c("sample" = "Col_ID")) |>
      filter(Group_Time %in% GROUP_ORDER) |>
      mutate(Group_Time = factor(Group_Time, levels = GROUP_ORDER))

    me_summary <- me_df |>
      group_by(Group_Time) |>
      summarise(mean_me = mean(me), se = sd(me) / sqrt(n()), .groups = "drop")

    p_me <- ggplot(me_summary, aes(x = Group_Time, y = mean_me, group = 1)) +
      geom_line(linewidth = 0.8, color = "grey30") +
      geom_pointrange(aes(ymin = mean_me - se, ymax = mean_me + se,
                          fill = Group_Time),
                      shape = 21, size = 1.5, color = "black") +
      scale_fill_manual(values = GROUP_COLS) +
      labs(title = "Eigengene", x = NULL, y = "ME") +
      FIG_THEME +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "none", plot.title = element_text(size = 9))
  } else {
    p_me <- ggplot() + theme_void() + labs(title = "ME not found")
  }

  # Right: ORA bars
  mod_enrich <- enrich_df |>
    filter(module == mod_color, padj < 0.05) |>
    arrange(padj) |>
    head(8) |>
    mutate(pathway_clean = clean_pathway_name(pathway))

  if (nrow(mod_enrich) > 0) {
    p_ora <- ggplot(mod_enrich,
                    aes(x = -log10(padj),
                        y = reorder(pathway_clean, -log10(padj)))) +
      geom_col(fill = "steelblue", width = 0.7) +
      labs(title = "Enrichment", x = "-log10(padj)", y = NULL) +
      FIG_THEME +
      theme(axis.text.y = element_text(size = 6),
            plot.title = element_text(size = 9))
  } else {
    p_ora <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "No sig pathways", size = 3)
  }

  p_heat | p_me | p_ora + plot_layout(widths = c(0.40, 0.25, 0.35))
}

# --- Build composite ---
rows <- lapply(key_mods, make_triptych_row)
pB <- wrap_plots(rows, ncol = 1)

PB_W <- 350; PB_H <- length(key_mods) * 120
ggsave(file.path(RPT, "panel_B_triptych_MAIN.pdf"), pB,
       width = PB_W, height = PB_H, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "panel_B_triptych_MAIN.png"), pB,
       width = PB_W, height = PB_H, units = "mm", dpi = 300,
       limitsize = FALSE)
message("F06 Panel B done")
