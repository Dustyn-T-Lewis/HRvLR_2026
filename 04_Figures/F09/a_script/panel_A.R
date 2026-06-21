# F09 Panel A: Module-Trait Correlation Heatmap (HR/LR-Stratified)
# Layout: gene counts (left) | 5-section heatmap
# Sections: Contrasts (LMM) | Baseline HR | Baseline LR | Delta HR | Delta LR
# BH correction: per-section for LMM, per-column for stratified
# Two-tier display: solid border = FDR < 0.05; dot = nominal p < 0.05
# Generates: panel_A_heatmap_MAIN.pdf/png, c_data/01_panel_A_*.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F09/a_script/style.R")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(patchwork); library(ggplot2)
})

RPT <- "04_Figures/F09/b_reports"

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

# --- Data loading ---

# LMM contrasts (lmer eigengene ~ group_time + (1|subject), KR df)
lmm_audit <- read_csv(file.path(WGCNA_DATA, "wgcna_lmm_contrast_audit.csv"),
                       show_col_types = FALSE)
lmm_r <- lmm_audit %>%
  select(module, contrast, r_equiv) %>%
  pivot_wider(names_from = contrast, values_from = r_equiv) %>%
  column_to_rownames("module") %>% as.matrix()
lmm_p <- lmm_audit %>%
  select(module, contrast, p_bh) %>%
  pivot_wider(names_from = contrast, values_from = p_bh) %>%
  column_to_rownames("module") %>% as.matrix()
lmm_p_raw <- lmm_audit %>%
  select(module, contrast, p_raw) %>%
  pivot_wider(names_from = contrast, values_from = p_raw) %>%
  column_to_rownames("module") %>% as.matrix()

# Stratified baseline correlations (HR and LR Pre separately, per-column BH)
bl_cor_hr <- read_csv(file.path(WGCNA_DATA, "wgcna_baseline_trait_correlations_hr.csv"),
                      show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
bl_pval_hr <- read_csv(file.path(WGCNA_DATA, "wgcna_baseline_trait_pvalues_bh_hr.csv"),
                       show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
bl_cor_lr <- read_csv(file.path(WGCNA_DATA, "wgcna_baseline_trait_correlations_lr.csv"),
                      show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
bl_pval_lr <- read_csv(file.path(WGCNA_DATA, "wgcna_baseline_trait_pvalues_bh_lr.csv"),
                       show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()

# Stratified change correlations (HR and LR paired separately, per-column BH)
ch_cor_hr <- read_csv(file.path(WGCNA_DATA, "wgcna_change_trait_correlations_hr.csv"),
                      show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
ch_pval_hr <- read_csv(file.path(WGCNA_DATA, "wgcna_change_trait_pvalues_bh_hr.csv"),
                       show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
ch_cor_lr <- read_csv(file.path(WGCNA_DATA, "wgcna_change_trait_correlations_lr.csv"),
                      show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
ch_pval_lr <- read_csv(file.path(WGCNA_DATA, "wgcna_change_trait_pvalues_bh_lr.csv"),
                       show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()

# Raw p-values for nominal significance tier
bl_praw_hr <- read_csv(file.path(WGCNA_DATA, "wgcna_baseline_trait_pvalues_raw_hr.csv"),
                       show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
bl_praw_lr <- read_csv(file.path(WGCNA_DATA, "wgcna_baseline_trait_pvalues_raw_lr.csv"),
                       show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
ch_praw_hr <- read_csv(file.path(WGCNA_DATA, "wgcna_change_trait_pvalues_raw_hr.csv"),
                       show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()
ch_praw_lr <- read_csv(file.path(WGCNA_DATA, "wgcna_change_trait_pvalues_raw_lr.csv"),
                       show_col_types = FALSE) %>%
  column_to_rownames("module") %>% as.matrix()

module_df      <- read_csv(file.path(WGCNA_DATA, "wgcna_module_assignments.csv"),
                           show_col_types = FALSE)
mod_bio_labels <- read_csv(file.path(PANEL_DATA, "mod_bio_labels.csv"),
                           show_col_types = FALSE)

if (!"display_label" %in% colnames(mod_bio_labels)) {
  mod_bio_labels <- mod_bio_labels %>%
    mutate(display_label = paste0(module_id, ": ", str_to_title(module_color)))
}
mod_display_vec <- setNames(mod_bio_labels$display_label, mod_bio_labels$module_color)
mod_display_label <- function(color) {
  lbl <- mod_display_vec[color]
  lbl[is.na(lbl)] <- str_to_title(color[is.na(lbl)])
  lbl
}

message("Panel A: five-section module-trait heatmap (HR/LR-stratified)...")

PA_W <- 420
PA_H <- 230

txt_cell  <- scale_text(BASE_GENE, PA_W) * 0.85
txt_count <- scale_text(BASE_COUNT, PA_W) * 0.85
txt_axis  <- scale_text(BASE_STAT, PA_W) * 1.5
txt_brack <- scale_text(BASE_GENE, PA_W) * 1.15

# --- Build combined 5-section correlation and p-value matrices (non-grey) ---
non_grey <- rownames(lmm_r)[rownames(lmm_r) != "MEgrey"]

lmm_cor  <- lmm_r[non_grey, , drop = FALSE]
lmm_pval <- lmm_p[non_grey, , drop = FALSE]

# Suffix stratified columns for uniqueness
bl_cor_h  <- bl_cor_hr[non_grey, , drop = FALSE]
bl_pval_h <- bl_pval_hr[non_grey, , drop = FALSE]
colnames(bl_cor_h)  <- paste0(colnames(bl_cor_h), "_HR")
colnames(bl_pval_h) <- paste0(colnames(bl_pval_h), "_HR")

bl_cor_l  <- bl_cor_lr[non_grey, , drop = FALSE]
bl_pval_l <- bl_pval_lr[non_grey, , drop = FALSE]
colnames(bl_cor_l)  <- paste0(colnames(bl_cor_l), "_LR")
colnames(bl_pval_l) <- paste0(colnames(bl_pval_l), "_LR")

ch_cor_h  <- ch_cor_hr[non_grey, , drop = FALSE]
ch_pval_h <- ch_pval_hr[non_grey, , drop = FALSE]
colnames(ch_cor_h)  <- paste0(colnames(ch_cor_h), "_HR")
colnames(ch_pval_h) <- paste0(colnames(ch_pval_h), "_HR")

ch_cor_l  <- ch_cor_lr[non_grey, , drop = FALSE]
ch_pval_l <- ch_pval_lr[non_grey, , drop = FALSE]
colnames(ch_cor_l)  <- paste0(colnames(ch_cor_l), "_LR")
colnames(ch_pval_l) <- paste0(colnames(ch_pval_l), "_LR")

cor_mat    <- cbind(lmm_cor, bl_cor_h, bl_cor_l, ch_cor_h, ch_cor_l)
pval_mat_b <- cbind(lmm_pval, bl_pval_h, bl_pval_l, ch_pval_h, ch_pval_l)

# Raw p-value matrix for nominal significance tier
lmm_pval_raw <- lmm_p_raw[non_grey, , drop = FALSE]
bl_praw_h <- bl_praw_hr[non_grey, , drop = FALSE]
colnames(bl_praw_h) <- paste0(colnames(bl_praw_h), "_HR")
bl_praw_l <- bl_praw_lr[non_grey, , drop = FALSE]
colnames(bl_praw_l) <- paste0(colnames(bl_praw_l), "_LR")
ch_praw_h <- ch_praw_hr[non_grey, , drop = FALSE]
colnames(ch_praw_h) <- paste0(colnames(ch_praw_h), "_HR")
ch_praw_l <- ch_praw_lr[non_grey, , drop = FALSE]
colnames(ch_praw_l) <- paste0(colnames(ch_praw_l), "_LR")
pval_raw_mat <- cbind(lmm_pval_raw, bl_praw_h, bl_praw_l, ch_praw_h, ch_praw_l)

trait_order <- colnames(cor_mat)
n_lmm <- ncol(lmm_cor)
n_bl_hr <- ncol(bl_cor_h)
n_bl_lr <- ncol(bl_cor_l)
n_ch_hr <- ncol(ch_cor_h)
n_ch_lr <- ncol(ch_cor_l)

# Clean column labels
col_labels <- setNames(
  c(
    gsub("_", "\n", colnames(lmm_cor)),
    gsub("_HR$", "", colnames(bl_cor_h)),
    gsub("_LR$", "", colnames(bl_cor_l)),
    gsub("_HR$", "", colnames(ch_cor_h)),
    gsub("_LR$", "", colnames(ch_cor_l))
  ),
  trait_order
)
# Shorten
col_labels <- gsub("X1RM_Leg_Pre", "1RM", col_labels)
col_labels <- gsub("mCSA_Pre", "mCSA", col_labels)
col_labels <- gsub("fCSA_Mixed_Pre", "fCSA_Mix", col_labels)
col_labels <- gsub("fCSA_Type_I_Pre", "fCSA_TI", col_labels)
col_labels <- gsub("fCSA_Type_II_Pre", "fCSA_TII", col_labels)
col_labels <- gsub("delta_X1RM", "\u03941RM", col_labels)
col_labels <- gsub("delta_mCSA", "\u0394mCSA", col_labels)
col_labels <- gsub("delta_fCSA_Mix", "\u0394fCSA_Mix", col_labels)
col_labels <- gsub("delta_fCSA_TI", "\u0394fCSA_TI", col_labels)
col_labels <- gsub("delta_fCSA_TII", "\u0394fCSA_TII", col_labels)
col_labels <- gsub("Baseline\nHRvLR", "Base.", col_labels)
col_labels <- gsub("Training\nHR", "Trn_HR", col_labels)
col_labels <- gsub("Training\nLR", "Trn_LR", col_labels)
col_labels <- gsub("Training\nInteraction", "Trn\u00d7Grp", col_labels)

# Column x-positions with gaps between sections
gap <- 0.3
n_total_cols <- n_lmm + n_bl_hr + n_bl_lr + n_ch_hr + n_ch_lr
xpos <- numeric(n_total_cols)
x <- 0
for (i in seq_len(n_lmm)) { x <- x + 1; xpos[i] <- x }
for (i in seq_len(n_bl_hr)) { x <- x + (if (i == 1) 1 + gap else 1); xpos[n_lmm + i] <- x }
for (i in seq_len(n_bl_lr)) { x <- x + (if (i == 1) 1 + gap else 1); xpos[n_lmm + n_bl_hr + i] <- x }
for (i in seq_len(n_ch_hr)) { x <- x + (if (i == 1) 1 + gap else 1); xpos[n_lmm + n_bl_hr + n_bl_lr + i] <- x }
for (i in seq_len(n_ch_lr)) { x <- x + (if (i == 1) 1 + gap else 1); xpos[n_lmm + n_bl_hr + n_bl_lr + n_ch_hr + i] <- x }
trait_xpos <- setNames(xpos, trait_order)

# Module order: hclust on correlation matrix
me_dist <- dist(cor_mat)
me_hc   <- hclust(me_dist)
mod_order <- rownames(cor_mat)[me_hc$order]

heat_df <- expand.grid(module = rownames(cor_mat), trait = colnames(cor_mat),
                       stringsAsFactors = FALSE) %>%
  mutate(cor      = as.vector(cor_mat),
         pval     = as.vector(pval_mat_b),
         pval_raw = as.vector(pval_raw_mat),
         stars    = sig_stars(pval),
         label    = sprintf("%.2f%s", cor, stars),
         trait_label = col_labels[trait],
         xpos     = trait_xpos[trait])
heat_df$module <- factor(heat_df$module, levels = mod_order)

mod_color_raw <- setNames(gsub("^ME", "", mod_order), mod_order)
mod_display   <- mod_display_label(mod_color_raw)
mod_display_map <- setNames(mod_display, mod_order)

light_modules <- c("yellow", "pink")

# --- Gene counts bar plot ---
mod_counts <- module_df %>%
  filter(module_color != "grey") %>%
  count(module_color, name = "n_proteins") %>%
  mutate(module = paste0("ME", module_color)) %>%
  filter(module %in% mod_order) %>%
  mutate(module = factor(module, levels = mod_order))

p_counts <- ggplot(mod_counts, aes(x = n_proteins, y = module)) +
  geom_col(fill = mod_counts$module_color, color = "black",
           linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = n_proteins, x = n_proteins / 2),
            size = txt_count, fontface = "bold",
            color = ifelse(mod_counts$module_color %in% light_modules,
                           "grey30", "white")) +
  scale_x_reverse(expand = c(0, 0), limits = c(max(mod_counts$n_proteins) * 1.02, 0)) +
  scale_y_discrete(labels = mod_display_map) +
  labs(y = NULL, x = "Gene Counts") +
  FIG_THEME +
  theme(axis.text.y        = element_text(size = txt_cell * 1.3, face = "bold"),
        axis.ticks.y       = element_blank(),
        axis.text.x        = element_text(size = txt_axis),
        axis.title.x       = element_text(size = txt_axis),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color = "black", linewidth = 0.3),
        legend.position    = "none",
        plot.margin        = margin(2, 0, 2, 2))

# --- Section brackets ---
xmin_all <- min(xpos) - 0.55
xmax_all <- max(xpos) + 0.55

# Compute shared objects for bracket building
shared_obj <- readRDS(file.path(PANEL_DATA, "shared_objects.rds"))
n_hr <- length(shared_obj$hr_subj)
n_lr <- length(shared_obj$lr_subj)
n_hr_paired <- length(shared_obj$common_hr)
n_lr_paired <- length(shared_obj$common_lr)

sec_lmm_start <- xpos[1]
sec_lmm_end   <- xpos[n_lmm]
sec_blhr_start <- xpos[n_lmm + 1]
sec_blhr_end   <- xpos[n_lmm + n_bl_hr]
sec_bllr_start <- xpos[n_lmm + n_bl_hr + 1]
sec_bllr_end   <- xpos[n_lmm + n_bl_hr + n_bl_lr]
sec_chhr_start <- xpos[n_lmm + n_bl_hr + n_bl_lr + 1]
sec_chhr_end   <- xpos[n_lmm + n_bl_hr + n_bl_lr + n_ch_hr]
sec_chlr_start <- xpos[n_lmm + n_bl_hr + n_bl_lr + n_ch_hr + 1]
sec_chlr_end   <- xpos[n_total_cols]

brackets <- tribble(
  ~label,                                                         ~start,         ~end,
  sprintf("Contrasts (LMM)\n(n=%d, paired)", n_hr + n_lr),       sec_lmm_start,  sec_lmm_end,
  sprintf("Baseline HR\n(n=%d Pre)", n_hr),                      sec_blhr_start, sec_blhr_end,
  sprintf("Baseline LR\n(n=%d Pre)", n_lr),                      sec_bllr_start, sec_bllr_end,
  sprintf("\u0394 HR\n(n=%d paired)", n_hr_paired),               sec_chhr_start, sec_chhr_end,
  sprintf("\u0394 LR\n(n=%d paired)", n_lr_paired),               sec_chlr_start, sec_chlr_end
)
brackets <- brackets %>% mutate(mid = (start + end) / 2)

p_brackets <- ggplot() +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = end + 0.4, y = 0.28, yend = 0.28),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = start - 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = end + 0.4, xend = end + 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  geom_text(data = brackets,
            aes(x = mid, y = 0.62, label = label),
            size = txt_brack, fontface = "bold", color = "grey25",
            lineheight = 0.85) +
  scale_x_continuous(limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 0, 0))

# --- Heatmap ---
p_heat <- ggplot(heat_df, aes(x = xpos, y = module, fill = cor)) +
  # Section background bands (group-tinted shading)
  annotate("rect", xmin = sec_lmm_start - 0.5, xmax = sec_lmm_end + 0.5,
           ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.06) +
  annotate("rect", xmin = sec_blhr_start - 0.5, xmax = sec_blhr_end + 0.5,
           ymin = -Inf, ymax = Inf, fill = GROUP_COLORS["HR"], alpha = 0.06) +
  annotate("rect", xmin = sec_bllr_start - 0.5, xmax = sec_bllr_end + 0.5,
           ymin = -Inf, ymax = Inf, fill = GROUP_COLORS["LR"], alpha = 0.06) +
  annotate("rect", xmin = sec_chhr_start - 0.5, xmax = sec_chhr_end + 0.5,
           ymin = -Inf, ymax = Inf, fill = GROUP_COLORS["HR"], alpha = 0.06) +
  annotate("rect", xmin = sec_chlr_start - 0.5, xmax = sec_chlr_end + 0.5,
           ymin = -Inf, ymax = Inf, fill = GROUP_COLORS["LR"], alpha = 0.06) +
  geom_tile(color = "black", linewidth = 0.3) +
  # Red/blue border for highest/lowest r per row
  geom_tile(data = heat_df %>% group_by(module) %>% filter(cor == max(cor)) %>% slice(1) %>% ungroup(),
            color = "#B2182B", linewidth = 1.3, fill = NA) +
  geom_tile(data = heat_df %>% group_by(module) %>% filter(cor == min(cor)) %>% slice(1) %>% ungroup(),
            color = "#2166AC", linewidth = 1.3, fill = NA) +
  # Solid thick black border for BH FDR < 0.05
  geom_tile(data = heat_df %>% filter(pval < 0.05),
            color = "black", linewidth = 1.0, fill = NA) +
  # Dot for nominal p < 0.05 (but FDR >= 0.05)
  geom_point(data = heat_df %>% filter(pval_raw < 0.05 & pval >= 0.05),
             shape = 16, size = 0.8, color = "grey30") +
  # Text labels
  geom_text(data = heat_df %>% filter(pval >= 0.05),
            aes(label = label), size = txt_cell, color = "black") +
  geom_text(data = heat_df %>% filter(pval < 0.05),
            aes(label = label), size = txt_cell * 1.15,
            fontface = "bold", color = "white") +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-0.8, 0.8),
                       oob = scales::squish,
                       name = "r / r-equiv.",
                       na.value = "white",
                       guide = guide_colorbar(
                         barwidth = unit(60, "mm"),
                         barheight = unit(4, "mm"),
                         title.position = "left",
                         title.vjust = 0.75
                       )) +
  scale_x_continuous(
    breaks = xpos,
    labels = col_labels[trait_order],
    limits = c(xmin_all, xmax_all),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1,
                                           size = txt_cell * 1.3, face = "bold"),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.ticks.x       = element_blank(),
        panel.grid         = element_blank(),
        panel.border       = element_blank(),
        legend.position    = "bottom",
        legend.title       = element_text(size = txt_axis, face = "bold"),
        legend.text        = element_text(size = txt_cell),
        plot.margin        = margin(0, 2, 2, -3))

# --- Layout assembly ---
sft_csv <- read_csv(file.path(WGCNA_DATA, "wgcna_sft_summary.csv"), show_col_types = FALSE)

design_layout <- c(
  area(1, 4, 1, 14),   # brackets (above heatmap)
  area(2, 1, 10, 3),   # gene counts
  area(2, 4, 10, 14)   # heatmap
)
pA <- wrap_elements(p_brackets) + p_counts + p_heat +
  plot_layout(design = design_layout)

pA <- pA +
  plot_annotation(
    title    = "WGCNA Module-Trait Associations",
    subtitle = sprintf("%d modules (Pearson, power %d) | LMM: lmer + KR df, BH-corrected | Stratified: bicor, BH per-trait",
                        length(mod_order), sft_csv$selected_power[1]),
    caption  = paste0("Values: r or r-equiv. | BH per-column (LMM: ", nrow(lmm_audit),
                      " tests; stratified: per trait) | ",
                      "Solid border = FDR < .05 | Dot = nominal p < .05 | ",
                      "Red border = highest r per row | Blue border = lowest r per row"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = txt_axis * 2.5),
      plot.subtitle = element_text(size = txt_cell * 2, color = "grey30", face = "italic"),
      plot.caption  = element_text(size = txt_cell * 1.5, color = "grey40", hjust = 0),
      plot.margin   = margin(2, 2, 2, 2)
    )
  )

write_csv(heat_df, file.path(DAT, "01_panel_A_heatmap_data.csv"))

ggsave(file.path(RPT_PDF, "panel_A_heatmap_MAIN.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "panel_A_heatmap_MAIN.png"), pA,
       width = PA_W, height = PA_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("  Panel A saved")
