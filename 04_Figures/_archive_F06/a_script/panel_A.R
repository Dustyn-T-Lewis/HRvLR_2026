# F06 Panel A: Module-Trait Heatmap (Five-Section, Group-Stratified)
# Layout: gene counts | 5-section heatmap (LMM | BL HR | BL LR | Δ HR | Δ LR)
# Two-tier display: solid border = FDR < 0.05; dashed = nominal p < 0.05

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

# --- Data loading ---
lmm_audit <- read_csv(file.path(DAT, "wgcna_lmm_contrast_audit.csv"), show_col_types = FALSE)
lmm_r <- lmm_audit |>
  select(module, contrast, r_equiv) |>
  pivot_wider(names_from = contrast, values_from = r_equiv) |>
  column_to_rownames("module") |> as.matrix()
lmm_p <- lmm_audit |>
  select(module, contrast, p_bh) |>
  pivot_wider(names_from = contrast, values_from = p_bh) |>
  column_to_rownames("module") |> as.matrix()
lmm_p_raw <- lmm_audit |>
  select(module, contrast, p_raw) |>
  pivot_wider(names_from = contrast, values_from = p_raw) |>
  column_to_rownames("module") |> as.matrix()

# Stratified baseline/change correlations
read_strat <- function(prefix, grp) {
  read_csv(file.path(DAT, sprintf("wgcna_%s_%s.csv", prefix, grp)),
           show_col_types = FALSE) |>
    column_to_rownames("module") |> as.matrix()
}

bl_cor_hr  <- read_strat("baseline_trait_correlations", "hr")
bl_pval_hr <- read_strat("baseline_trait_pvalues_bh", "hr")
bl_cor_lr  <- read_strat("baseline_trait_correlations", "lr")
bl_pval_lr <- read_strat("baseline_trait_pvalues_bh", "lr")
ch_cor_hr  <- read_strat("change_trait_correlations", "hr")
ch_pval_hr <- read_strat("change_trait_pvalues_bh", "hr")
ch_cor_lr  <- read_strat("change_trait_correlations", "lr")
ch_pval_lr <- read_strat("change_trait_pvalues_bh", "lr")

bl_praw_hr <- read_strat("baseline_trait_pvalues_raw", "hr")
bl_praw_lr <- read_strat("baseline_trait_pvalues_raw", "lr")
ch_praw_hr <- read_strat("change_trait_pvalues_raw", "hr")
ch_praw_lr <- read_strat("change_trait_pvalues_raw", "lr")

# Module metadata
module_df <- read_csv(file.path(DAT, "wgcna_module_assignments.csv"), show_col_types = FALSE)
gene_counts <- module_df |>
  count(module_color, name = "n_genes") |>
  filter(module_color != "grey")

modules <- rownames(lmm_r)
modules <- modules[modules != "grey"]

# --- Build tile data for each section ---
build_section <- function(r_mat, p_mat, p_raw_mat, section_name) {
  mods <- intersect(modules, rownames(r_mat))
  r_mat <- r_mat[mods, , drop = FALSE]
  p_mat <- p_mat[mods, , drop = FALSE]
  p_raw_mat <- p_raw_mat[mods, , drop = FALSE]

  expand_grid(module = mods, trait = colnames(r_mat)) |>
    mutate(
      r     = mapply(function(m, t) r_mat[m, t], module, trait),
      p_bh  = mapply(function(m, t) p_mat[m, t], module, trait),
      p_raw = mapply(function(m, t) p_raw_mat[m, t], module, trait),
      sig   = case_when(p_bh < 0.05 ~ "FDR", p_raw < 0.05 ~ "nom", TRUE ~ "ns"),
      section = section_name
    )
}

tile_lmm <- build_section(lmm_r, lmm_p, lmm_p_raw, "LMM Contrasts")
tile_bl_hr <- build_section(bl_cor_hr, bl_pval_hr, bl_praw_hr, "Baseline HR")
tile_bl_lr <- build_section(bl_cor_lr, bl_pval_lr, bl_praw_lr, "Baseline LR")
tile_ch_hr <- build_section(ch_cor_hr, ch_pval_hr, ch_praw_hr, "\u0394 HR")
tile_ch_lr <- build_section(ch_cor_lr, ch_pval_lr, ch_praw_lr, "\u0394 LR")

tile_all <- bind_rows(tile_lmm, tile_bl_hr, tile_bl_lr, tile_ch_hr, tile_ch_lr) |>
  mutate(
    section = factor(section, levels = c("LMM Contrasts", "Baseline HR",
                                          "Baseline LR", "\u0394 HR", "\u0394 LR")),
    module  = factor(module, levels = rev(modules))
  )

max_abs <- max(abs(tile_all$r), na.rm = TRUE)

# --- Plot ---
p_heat <- ggplot(tile_all, aes(x = trait, y = module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(data = tile_all |> filter(sig == "FDR"),
            color = "black", linewidth = 0.8, fill = NA) +
  geom_tile(data = tile_all |> filter(sig == "nom"),
            color = "black", linewidth = 0.4, linetype = "dashed", fill = NA) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-max_abs, max_abs),
                       name = "r") +
  facet_grid(~ section, scales = "free_x", space = "free_x") +
  labs(title = "Module-Trait Associations",
       subtitle = "Solid border = FDR < 0.05 | Dashed = p < 0.05",
       x = NULL, y = NULL, tag = "A") +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 8, face = "bold"),
    strip.text  = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    legend.key.width = unit(15, "mm")
  )

# Gene count bar
gc_df <- gene_counts |>
  filter(module_color %in% modules) |>
  mutate(module_color = factor(module_color, levels = rev(modules)))

p_counts <- ggplot(gc_df, aes(x = n_genes, y = module_color)) +
  geom_col(fill = "grey60", width = 0.7) +
  geom_text(aes(label = n_genes), hjust = -0.1, size = 2.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "Genes", y = NULL) +
  FIG_THEME +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank())

pA <- p_heat | p_counts + plot_layout(widths = c(0.85, 0.15))

PA_W <- 350; PA_H <- max(150, length(modules) * 15 + 60)
ggsave(file.path(RPT, "panel_A_module_trait_MAIN.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_A_module_trait_MAIN.png"), pA,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
message("F06 Panel A done")
