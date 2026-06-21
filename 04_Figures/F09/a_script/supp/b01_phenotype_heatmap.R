# Supplementary: Module-Phenotype Heatmap (all modules, all traits, pooled)
# Complements Panel A (which is stratified by HR/LR). This shows pooled
# baseline + change correlations for all modules, unstratified.
# Generates: b01_phenotype_heatmap_SUPP.pdf/png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F09/a_script/style.R")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(ggplot2); library(patchwork); library(WGCNA)
})

RPT <- "04_Figures/F09/b_reports/supp"

RPT_PDF       <- file.path(RPT, "main", "pdf")

RPT_PNG       <- file.path(RPT, "main", "png")

RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")

RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
DAT <- "04_Figures/F09/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("Supplementary: module-phenotype heatmap (pooled)...")

# --- Load data ---
MEs        <- readRDS(file.path(PANEL_DATA, "MEs.rds"))
meta       <- read_csv(file.path(PANEL_DATA, "meta.csv"), show_col_types = FALSE)
me_pre     <- readRDS(file.path(PANEL_DATA, "me_pre.rds"))
me_post    <- readRDS(file.path(PANEL_DATA, "me_post.rds"))
delta_me   <- readRDS(file.path(PANEL_DATA, "delta_me.rds"))
pheno_wide <- read_csv(file.path(PANEL_DATA, "pheno_wide.csv"), show_col_types = FALSE)
subj_group <- read_csv(file.path(PANEL_DATA, "subj_group.csv"), show_col_types = FALSE)
module_df  <- read_csv(file.path(WGCNA_DATA, "wgcna_module_assignments.csv"),
                        show_col_types = FALSE)
mod_bio    <- read_csv(file.path(PANEL_DATA, "mod_bio_labels.csv"), show_col_types = FALSE)

if (!"display_label" %in% colnames(mod_bio)) {
  mod_bio <- mod_bio %>%
    mutate(display_label = paste0(module_id, ": ", str_to_title(module_color)))
}
mod_display_vec <- setNames(mod_bio$display_label, mod_bio$module_color)

# --- Baseline correlations (pooled T1) ---
t1_meta <- meta %>% filter(timepoint == "T1")
t1_ids  <- t1_meta$sample_id
me_t1   <- MEs[t1_ids, , drop = FALSE]

baseline_traits <- data.frame(
  X1RM_Leg   = t1_meta$X1RM_Leg_Pre,
  mCSA       = t1_meta$mCSA_Pre,
  fCSA_Mixed = t1_meta$fCSA_Mixed_Pre,
  fCSA_TI    = t1_meta$fCSA_Type_I_Pre,
  fCSA_TII   = t1_meta$fCSA_Type_II_Pre,
  row.names  = t1_ids
)

all_mods <- colnames(MEs)
non_grey <- all_mods[all_mods != "MEgrey"]

bl_cor  <- bicor(me_t1[, non_grey], baseline_traits, use = "pairwise.complete.obs")
bl_pval <- matrix(NA_real_, nrow = length(non_grey), ncol = ncol(baseline_traits),
                  dimnames = list(non_grey, colnames(baseline_traits)))
for (trait in colnames(baseline_traits)) {
  n_ok <- sum(!is.na(baseline_traits[, trait]))
  bl_pval[, trait] <- corPvalueStudent(bl_cor[, trait], n_ok)
}

# --- Change correlations (pooled delta) ---
common_subj <- intersect(rownames(delta_me), pheno_wide$subject)
dme_pooled  <- delta_me[common_subj, non_grey, drop = FALSE]

change_traits <- data.frame(
  delta_X1RM     = pheno_wide$delta_X1RM[match(common_subj, pheno_wide$subject)],
  delta_mCSA     = pheno_wide$delta_mCSA[match(common_subj, pheno_wide$subject)],
  delta_fCSA_Mix = pheno_wide$delta_fCSA_Mix[match(common_subj, pheno_wide$subject)],
  delta_fCSA_TI  = pheno_wide$delta_fCSA_TI[match(common_subj, pheno_wide$subject)],
  delta_fCSA_TII = pheno_wide$delta_fCSA_TII[match(common_subj, pheno_wide$subject)],
  row.names = common_subj
)

ch_cor  <- bicor(dme_pooled, change_traits, use = "pairwise.complete.obs")
ch_pval <- matrix(NA_real_, nrow = length(non_grey), ncol = ncol(change_traits),
                  dimnames = list(non_grey, colnames(change_traits)))
for (trait in colnames(change_traits)) {
  n_ok <- sum(!is.na(change_traits[, trait]))
  ch_pval[, trait] <- corPvalueStudent(ch_cor[, trait], n_ok)
}

# --- Combined matrix ---
# BH per-column
bh_col <- function(pmat) {
  out <- pmat
  for (j in seq_len(ncol(pmat))) out[, j] <- p.adjust(pmat[, j], method = "BH")
  out
}
bl_pval_bh <- bh_col(bl_pval)
ch_pval_bh <- bh_col(ch_pval)

cor_all  <- cbind(bl_cor, ch_cor)
pval_all <- cbind(bl_pval_bh, ch_pval_bh)

# Module order by hclust
me_dist <- dist(cor_all)
me_hc   <- hclust(me_dist)
mod_order <- rownames(cor_all)[me_hc$order]

# Trait labels
trait_labels <- c(
  X1RM_Leg   = "1RM Leg",
  mCSA       = "mCSA",
  fCSA_Mixed = "fCSA Mixed",
  fCSA_TI    = "fCSA Type I",
  fCSA_TII   = "fCSA Type II",
  delta_X1RM     = "\u03941RM",
  delta_mCSA     = "\u0394mCSA",
  delta_fCSA_Mix = "\u0394fCSA Mix",
  delta_fCSA_TI  = "\u0394fCSA TI",
  delta_fCSA_TII = "\u0394fCSA TII"
)

heat_df <- expand.grid(module = rownames(cor_all), trait = colnames(cor_all),
                       stringsAsFactors = FALSE) %>%
  mutate(cor    = as.vector(cor_all),
         pval   = as.vector(pval_all),
         stars  = sig_stars(pval),
         label  = sprintf("%.2f%s", cor, stars),
         module = factor(module, levels = mod_order),
         trait_label = trait_labels[trait],
         trait  = factor(trait, levels = colnames(cor_all)))

mod_color_raw <- setNames(gsub("^ME", "", mod_order), mod_order)
mod_display   <- setNames(
  ifelse(mod_color_raw %in% names(mod_display_vec),
         mod_display_vec[mod_color_raw],
         str_to_title(mod_color_raw)),
  mod_order
)

PH_W <- 280
PH_H <- 200
txt_cell <- scale_text(BASE_GENE, PH_W) * 0.85
txt_axis <- scale_text(BASE_STAT, PH_W) * 1.3

p_pheno <- ggplot(heat_df, aes(x = trait, y = module, fill = cor)) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_tile(data = heat_df %>% filter(pval < 0.05),
            color = "black", linewidth = 1.0, fill = NA) +
  geom_text(data = heat_df %>% filter(pval >= 0.05),
            aes(label = label), size = txt_cell, color = "black") +
  geom_text(data = heat_df %>% filter(pval < 0.05),
            aes(label = label), size = txt_cell * 1.15,
            fontface = "bold", color = "white") +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-0.8, 0.8),
                       oob = scales::squish,
                       name = "bicor",
                       guide = guide_colorbar(barwidth = unit(50, "mm"),
                                              barheight = unit(4, "mm"))) +
  scale_x_discrete(labels = trait_labels[colnames(cor_all)]) +
  scale_y_discrete(labels = mod_display) +
  labs(title = "Module-Phenotype Correlations (Pooled, All Modules)",
       subtitle = sprintf("n=%d baseline | n=%d paired | bicor, BH per-trait | Solid border = FDR < .05",
                           nrow(t1_meta), length(common_subj)),
       x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = txt_axis, face = "bold"),
        axis.text.y = element_text(size = txt_cell * 1.2, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = txt_axis, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_blank())

ggsave(file.path(RPT_SUPP_PDF, "b01_phenotype_heatmap_SUPP.pdf"), p_pheno,
       width = PH_W, height = PH_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_SUPP_PNG, "b01_phenotype_heatmap_SUPP.png"), p_pheno,
       width = PH_W, height = PH_H, units = "mm",
       dpi = 300, limitsize = FALSE)

write_csv(heat_df, file.path(DAT, "b01_phenotype_heatmap_data.csv"))

message("  Supplementary phenotype heatmap saved")
