# F07 Panel D: fGSEA NES Scatter (GO Slim + Hallmark, a priori) -- Concordance
# NES scatter: Training_HR NES vs Training_LR NES
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F07/a_script/style.R")

library(tidyverse)
library(ggrepel)

PG_W <- 200
RPT <- "04_Figures/F07/b_reports"
DAT <- "04_Figures/F07/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# -- Load fGSEA cache from F05 (DEP figure that generated it) -----------------
fgsea_cache <- "04_Figures/F05/c_data/fgsea_cache.csv"
if (!file.exists(fgsea_cache))
  stop("fGSEA cache not found -- run F05 panel scripts first")

fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

# Filter to a priori collection: GO:BP (as GO Slim proxy) + Hallmark
fgsea_hg <- fgsea_all %>%
  filter(database %in% c("Hallmark", "GO:BP"),
         contrast %in% c("Training_HR", "Training_LR"))

# Pivot ALL terms to wide format
fgsea_wide <- fgsea_hg %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Training_HR), !is.na(NES_Training_LR)) %>%
  mutate(set_size = coalesce(size_Training_HR, size_Training_LR))

# Classify significance
fgsea_wide <- fgsea_wide %>%
  mutate(
    sig_HR = !is.na(padj_Training_HR) & padj_Training_HR < 0.05,
    sig_LR = !is.na(padj_Training_LR) & padj_Training_LR < 0.05,
    significance = case_when(
      sig_HR & sig_LR ~ "Sig Both",
      sig_HR          ~ "Sig HR only",
      sig_LR          ~ "Sig LR only",
      TRUE            ~ "NS"
    ) %>% factor(levels = names(SIG_COLORS_F7)),
    pathway_label = clean_pathway_name(pathway),
    db_shape = ifelse(database == "Hallmark", 24, 21)
  )

fgsea_sig <- fgsea_wide %>% filter(significance != "NS")

message(sprintf("  %d total pathways (Hallmark: %d, GO:BP: %d) | %d significant",
                nrow(fgsea_wide),
                sum(fgsea_wide$database == "Hallmark"),
                sum(fgsea_wide$database == "GO:BP"),
                nrow(fgsea_sig)))

# Spearman correlation
nes_cor_all <- cor.test(fgsea_wide$NES_Training_HR, fgsea_wide$NES_Training_LR, method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig$NES_Training_HR, fgsea_sig$NES_Training_LR, method = "spearman")
} else NULL

nes_lim <- max(abs(c(fgsea_wide$NES_Training_HR, fgsea_wide$NES_Training_LR))) * 1.15

# Quadrant counts (sig terms only)
n_conc_tr  <- sum(fgsea_sig$NES_Training_HR > 0 & fgsea_sig$NES_Training_LR > 0)
n_conc_bl  <- sum(fgsea_sig$NES_Training_HR < 0 & fgsea_sig$NES_Training_LR < 0)
n_disc_q2  <- sum(fgsea_sig$NES_Training_HR < 0 & fgsea_sig$NES_Training_LR > 0)
n_disc_q4  <- sum(fgsea_sig$NES_Training_HR > 0 & fgsea_sig$NES_Training_LR < 0)
n_conc_pw  <- n_conc_tr + n_conc_bl
n_total_sig <- nrow(fgsea_sig)
pw_conc_frac <- if (n_total_sig > 0) n_conc_pw / n_total_sig else 0

message(sprintf("  NES Spearman (all): rho = %.3f [%.3f, %.3f]",
                nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2]))

# Sizes
txt_pw   <- scale_text(BASE_PATHWAY, PG_W) * 1.25
txt_quad <- scale_text(BASE_QUADRANT, PG_W) * 1.4

# Label data for all sig terms
label_pw <- fgsea_sig %>%
  mutate(
    label_fill     = SIG_LABEL_FILL_F7[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F7[as.character(significance)],
    pathway_label  = pathway_label %>%
      str_replace("Amino Acid Metabolic.*", "Amino Acid Metabolism") %>%
      str_replace("Muscle System.*", "Muscle System") %>%
      str_replace("Ketone Metabolic.*", "Ketone Metabolism")
  )

# Split for layered plotting
ns_df  <- fgsea_wide %>% filter(significance == "NS")
sig_df <- fgsea_wide %>% filter(significance != "NS") %>%
  mutate(draw_order = factor(significance,
    levels = c("Sig LR only", "Sig HR only", "Sig Both", "Interaction"))) %>%
  arrange(draw_order)

# Subtitle
rho_sig_str <- if (!is.null(nes_cor_sig)) sprintf(", \u03c1(sig) = %.2f", nes_cor_sig$estimate) else ""
subtitle_str <- sprintf(
  "GO:BP + Hallmark (a priori) | %d pathways (%d sig.) | fGSEA on limma t-statistics\n\u03c1(all) = %.2f [%.2f, %.2f], %s%s | %.0f%% concordant",
  nrow(fgsea_wide), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", nes_cor_all$p.value)),
  rho_sig_str, pw_conc_frac * 100
)

pG <- ggplot(mapping = aes(x = NES_Training_HR, y = NES_Training_LR)) +
  # Quadrant backgrounds
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # NS points
  geom_point(data = ns_df, aes(shape = database),
             size = 2.5, fill = "grey70", color = "grey55", alpha = 0.40, stroke = 0.4) +
  # Sig points
  geom_point(data = sig_df, aes(fill = significance, size = set_size, shape = database),
             color = ifelse(sig_df$database == "Hallmark", "black", "grey65"),
             alpha = 0.80, stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS_F7, name = "Significance") +
  scale_shape_manual(values = c("Hallmark" = 24, "GO:BP" = 21), name = "Database") +
  scale_size_continuous(range = c(3, 10), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  # Pathway labels
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   size = txt_pw, fontface = "bold",
                   max.overlaps = 50,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, point.padding = 0.4,
                   force = 5, force_pull = 0.3,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  # Quadrant labels
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up  n = %d", n_conc_tr),
           hjust = 1.05, vjust = 1.2, size = txt_quad, fontface = "bold",
           color = "#D6604D", fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down  n = %d", n_conc_bl),
           hjust = -0.05, vjust = -0.2, size = txt_quad, fontface = "bold",
           color = "#D6604D", fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant  n = %d", n_disc_q2),
           hjust = -0.05, vjust = 1.2, size = txt_quad, fontface = "bold",
           color = "#4393C3", fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant  n = %d", n_disc_q4),
           hjust = 1.05, vjust = -0.2, size = txt_quad, fontface = "bold",
           color = "#4393C3", fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(title = "Pathway-Level Concordance (fGSEA)",
       subtitle = subtitle_str,
       x = "NES (Training HR)",
       y = "NES (Training LR)") +
  FIG_THEME +
  theme(
    axis.text      = element_text(size = 10, face = "bold", color = "grey30"),
    axis.title     = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(4, "mm"),
    legend.margin   = margin(0, 0, 0, 0),
    legend.box      = "horizontal"
  ) +
  guides(fill  = guide_legend(nrow = 1, order = 1,
                               override.aes = list(size = 4, alpha = 0.8, shape = 21)),
         shape = guide_legend(nrow = 1, order = 2,
                               override.aes = list(size = 4, fill = "grey50")),
         size  = guide_legend(nrow = 1, order = 3))

ggsave(file.path(RPT, "panel_D_nes_scatter_MAIN.pdf"), pG,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_nes_scatter_MAIN.png"), pG,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)

# Export all terms
fgsea_wide %>%
  transmute(
    pathway, pathway_label, database,
    NES_Training_HR = round(NES_Training_HR, 3),
    NES_Training_LR = round(NES_Training_LR, 3),
    padj_Training_HR = signif(padj_Training_HR, 4),
    padj_Training_LR = signif(padj_Training_LR, 4),
    significance     = as.character(significance),
    set_size
  ) %>%
  arrange(significance, desc(abs(NES_Training_HR) + abs(NES_Training_LR))) %>%
  write_csv(file.path(DAT, "panel_D", "nes_scatter.csv"))

cat("F07 Panel D done\n")
