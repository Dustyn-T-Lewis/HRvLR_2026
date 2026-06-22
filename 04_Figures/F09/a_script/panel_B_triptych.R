# F09 Panel B: WGCNA Per-Module Triptych
# Layout per row: z-score heatmap | eigengene dynamics | ORA bars
# Key modules ordered by LMM significance (from key_modules.txt)
# Eigengene dynamics: paired T1->T2 by HR/LR with LMM brackets

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F09/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

RPT <- "04_Figures/F09/b_reports"

RPT_PDF <- file.path(RPT, "main", "pdf")

RPT_PNG <- file.path(RPT, "main", "png")

RPT_SUPP_PDF <- file.path(RPT, "supp", "pdf")

RPT_SUPP_PNG <- file.path(RPT, "supp", "png")
DAT <- "04_Figures/F09/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("Panel B: WGCNA per-module triptych...")

# --- Load data ---
MEs <- readRDS(file.path(PANEL_DATA, "MEs.rds"))
datExpr <- readRDS(file.path(PANEL_DATA, "datExpr.rds"))
meta <- read_csv(file.path(PANEL_DATA, "meta.csv"), show_col_types = FALSE)
module_df <- read_csv(file.path(PANEL_DATA, "module_df.csv"), show_col_types = FALSE)
mod_bio <- read_csv(file.path(PANEL_DATA, "mod_bio_labels.csv"), show_col_types = FALSE)
lmm_audit <- read_csv(file.path(WGCNA_DATA, "wgcna_lmm_contrast_audit.csv"),
  show_col_types = FALSE
)
enrich_df <- read_csv(file.path(WGCNA_DATA, "wgcna_module_enrichment.csv"),
  show_col_types = FALSE
)
ann <- read_csv(file.path(PANEL_DATA, "imp_annotations.csv"), show_col_types = FALSE)

if (!"display_label" %in% colnames(mod_bio)) {
  mod_bio <- mod_bio %>%
    mutate(display_label = paste0(module_id, ": ", str_to_title(module_color)))
}
mod_labels <- setNames(mod_bio$display_label, mod_bio$module_color)

# Key modules
km_file <- file.path(PANEL_DATA, "key_modules.txt")
KEY_MODULES <- if (file.exists(km_file)) {
  readLines(km_file)
} else {
  km_file2 <- file.path(WGCNA_DATA, "key_modules.txt")
  if (file.exists(km_file2)) readLines(km_file2) else stop("key_modules.txt not found")
}
KEY_MODULES <- KEY_MODULES[nzchar(trimws(KEY_MODULES))]
mod_order <- KEY_MODULES

# --- Prepare z-score and eigengene data ---

# Z-score: per-protein z across samples, grouped by group_time
imp_mat <- readRDS(file.path(PANEL_DATA, "imp_mat.rds"))

group_order <- c("HR_T1", "HR_T2", "LR_T1", "LR_T2")
group_labels <- c(
  HR_T1 = "HR-T1", HR_T2 = "HR-T2",
  LR_T1 = "LR-T1", LR_T2 = "LR-T2"
)

# Build z-scores for each key module
z_list <- list()
for (mod in KEY_MODULES) {
  mod_ids <- module_df$uniprot_id[module_df$module_color == mod]
  mod_ids <- intersect(mod_ids, rownames(imp_mat))
  if (length(mod_ids) == 0) next

  sub_mat <- imp_mat[mod_ids, meta$sample_id, drop = FALSE]
  z_mat <- t(scale(t(sub_mat))) # z-score per protein across samples

  for (gt in group_order) {
    samps <- meta$sample_id[meta$group_time == gt]
    samps <- intersect(samps, colnames(z_mat))
    if (length(samps) == 0) next
    mean_z <- rowMeans(z_mat[, samps, drop = FALSE], na.rm = TRUE)
    genes <- ann$gene[match(names(mean_z), ann$uniprot_id)]
    z_list <- c(z_list, list(tibble(
      module = mod, uniprot_id = names(mean_z), gene = genes,
      group = gt, z = mean_z
    )))
  }
}
z_scores <- bind_rows(z_list)

# Build eigengene data
me_list <- list()
for (i in seq_len(nrow(meta))) {
  sid <- meta$sample_id[i]
  if (!(sid %in% rownames(MEs))) next
  for (mod in KEY_MODULES) {
    me_col <- paste0("ME", mod)
    if (!(me_col %in% colnames(MEs))) next
    me_list <- c(me_list, list(tibble(
      module     = mod,
      sample_id  = sid,
      subject    = meta$subject[i],
      group      = meta$group[i],
      timepoint  = meta$timepoint[i],
      group_time = meta$group_time[i],
      eigengene  = MEs[sid, me_col]
    )))
  }
}
me_data <- bind_rows(me_list)

# Save triptych intermediates
write_csv(z_scores, file.path(DAT, "03_panel_B_heatmap_zscores.csv"))
write_csv(me_data, file.path(DAT, "03_panel_B_eigengene_data.csv"))

# --- Dimensions ---
PB_W <- 280
PB_H <- 80 * length(mod_order) + 30

txt_heat <- scale_text(BASE_GENE, PB_W) * 0.7
txt_axis <- scale_text(BASE_STAT, PB_W) * 1.0
txt_title <- scale_text(BASE_GENE, PB_W) * 1.3
txt_bar <- scale_text(BASE_GENE, PB_W) * 0.95
txt_sig <- scale_text(BASE_GENE, PB_W) * 0.85

light_modules <- c("yellow", "pink")

# LMM stats for bracket annotations
lmm_stats <- lmm_audit %>%
  filter(contrast %in% c("Training_HR", "Training_LR", "Baseline_HRvLR")) %>%
  mutate(module = gsub("^ME", "", module)) %>%
  select(module, contrast, p_bh)

# --- Build one triptych row per module ---
build_row <- function(mod, show_xlab = FALSE) {
  label <- if (mod %in% names(mod_labels)) mod_labels[mod] else str_to_title(mod)
  n_mod <- module_df %>%
    filter(module_color == mod) %>%
    nrow()
  title_txt <- paste0(label, " (n=", n_mod, ")")

  # -- Heatmap: z-scores (4 group columns) --
  z_mod <- z_scores %>%
    filter(module == mod) %>%
    mutate(group = factor(group, levels = group_order))

  gene_order <- z_mod %>%
    filter(group == "HR_T1") %>%
    arrange(z) %>%
    pull(gene)
  z_mod$gene <- factor(z_mod$gene, levels = gene_order)

  p_heat <- ggplot(z_mod, aes(x = group, y = gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-2, 2), oob = scales::squish,
      guide = "none"
    ) +
    scale_x_discrete(labels = group_labels, position = "bottom") +
    labs(title = title_txt, y = NULL, x = NULL) +
    FIG_THEME +
    theme(
      plot.title = element_text(size = txt_title, face = "bold"),
      axis.text.x = if (show_xlab) {
        element_text(size = txt_axis * 0.85, angle = 45, hjust = 1)
      } else {
        element_blank()
      },
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(2, 1, 2, 2)
    )

  # -- Eigengene dynamics: paired T1->T2 with LMM brackets --
  me_mod <- me_data %>%
    filter(module == mod) %>%
    mutate(
      timepoint = factor(timepoint, levels = c("T1", "T2")),
      group = factor(group, levels = c("HR", "LR"))
    )

  # LMM p-values for brackets
  p_thr <- lmm_stats %>%
    filter(module == mod, contrast == "Training_HR") %>%
    pull(p_bh)
  p_tlr <- lmm_stats %>%
    filter(module == mod, contrast == "Training_LR") %>%
    pull(p_bh)
  p_bl <- lmm_stats %>%
    filter(module == mod, contrast == "Baseline_HRvLR") %>%
    pull(p_bh)

  fmt_sig <- function(p) {
    if (length(p) == 0 || is.na(p)) {
      return("ns")
    }
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  }

  p_eigen <- ggplot(me_mod, aes(x = timepoint, y = eigengene)) +
    geom_line(aes(group = subject, color = group), alpha = 0.2, linewidth = 0.3) +
    stat_summary(aes(group = group, color = group), fun = mean, geom = "line", linewidth = 1.2) +
    stat_summary(aes(group = group, color = group), fun = mean, geom = "point", size = 2.5) +
    scale_color_manual(values = GROUP_COLORS, guide = "none") +
    # Bracket annotations
    annotate("text",
      x = 1.5, y = max(me_mod$eigengene, na.rm = TRUE) * 0.95,
      label = fmt_sig(p_thr), size = txt_sig, fontface = "bold",
      color = GROUP_COLORS["HR"]
    ) +
    annotate("text",
      x = 1.5, y = min(me_mod$eigengene, na.rm = TRUE) * 0.95,
      label = fmt_sig(p_tlr), size = txt_sig, fontface = "bold",
      color = GROUP_COLORS["LR"]
    ) +
    annotate("text",
      x = 0.65,
      y = mean(c(
        max(me_mod$eigengene, na.rm = TRUE),
        min(me_mod$eigengene, na.rm = TRUE)
      )),
      label = fmt_sig(p_bl), size = txt_sig, fontface = "bold",
      color = "grey40", angle = 90
    ) +
    labs(y = "Eigengene", x = NULL) +
    FIG_THEME +
    theme(
      axis.text.x = if (show_xlab) element_text(size = txt_axis * 0.85) else element_blank(),
      axis.text.y = element_text(size = txt_axis * 0.75),
      axis.title.y = element_text(size = txt_axis * 0.8),
      panel.border = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.margin = margin(2, 1, 2, 1)
    )

  # -- ORA bars: top 5 per module --
  bar_data <- enrich_df %>%
    filter(module == mod, p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(5) %>%
    mutate(
      neg_log10_p = -log10(p.adjust),
      clean_name  = clean_pathway_name(Description),
      db_fill     = if ("database" %in% names(.)) DB_COLORS[database] else "#AA336A"
    ) %>%
    mutate(clean_name = factor(clean_name, levels = rev(clean_name)))

  if (nrow(bar_data) == 0) {
    p_bars <- ggplot() +
      annotate("text",
        x = 0.5, y = 0.5, label = "No sig.\nenrichment",
        size = txt_bar, color = "grey50"
      ) +
      theme_void() +
      theme(plot.margin = margin(2, 2, 2, 1))
  } else {
    p_bars <- ggplot(bar_data, aes(x = neg_log10_p, y = clean_name)) +
      geom_col(aes(fill = db_fill), color = "black", linewidth = 0.3, width = 0.7) +
      geom_text(aes(label = clean_name, x = 0.3),
        hjust = 0, size = txt_bar,
        fontface = "bold",
        color = ifelse(bar_data$db_fill %in% c("#AA336A", "#1565C0", "#00796B"),
          "white", "grey20"
        )
      ) +
      scale_fill_identity() +
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.05)),
        breaks = scales::breaks_pretty(n = 3),
        name = if (show_xlab) expression(-log[10](p[adj])) else NULL
      ) +
      scale_y_discrete(labels = NULL) +
      labs(y = NULL) +
      FIG_THEME +
      theme(
        axis.text.x  = if (show_xlab) element_text(size = txt_axis * 0.75) else element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x  = element_line(color = "black", linewidth = 0.3),
        panel.grid   = element_blank(),
        plot.margin  = margin(2, 2, 2, 1)
      )
  }

  # Combine row
  p_heat + p_eigen + p_bars + plot_layout(widths = c(3, 2, 3))
}

# --- Build all rows ---
rows <- lapply(seq_along(mod_order), function(i) {
  build_row(mod_order[i], show_xlab = (i == length(mod_order)))
})

# Z-score legend
z_legend <- ggplot(
  data.frame(z = seq(-2, 2, length.out = 100)),
  aes(x = z, y = 1, fill = z)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-2, 2),
    name = "Z-score", guide = guide_colorbar(
      barwidth = unit(40, "mm"), barheight = unit(3, "mm")
    )
  ) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = txt_axis * 0.7))

# --- Assemble ---
triptych <- wrap_plots(rows, ncol = 1) / wrap_elements(z_legend) +
  plot_layout(heights = c(rep(1, length(mod_order)), 0.12))

ggsave(file.path(RPT, "panel_b_triptych.pdf"), triptych,
  width = PB_W, height = PB_H, units = "mm",
  device = pdf_device, limitsize = FALSE, bg = "white"
)
ggsave(file.path(RPT, "panel_b_triptych.png"), triptych,
  width = PB_W, height = PB_H, units = "mm",
  dpi = 300, limitsize = FALSE, bg = "white"
)

message("  Panel B (WGCNA triptych) saved")
