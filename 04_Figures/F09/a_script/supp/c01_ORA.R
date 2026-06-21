# F09 Supp C01: ORA Scatter — Acute Concordance
# Acute concordance: Acute_HR (x) vs Acute_LR (y)
# (Moved from main panel to supplementary)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F09/a_script/style.R")
source("04_Figures/shared/concordance_utils.R")

suppressPackageStartupMessages({
  library(tidyverse); library(fgsea); library(ggrepel); library(patchwork)
})

RPT <- "04_Figures/F09/b_reports/supp"

RPT_PDF       <- file.path(RPT, "main", "pdf")

RPT_PNG       <- file.path(RPT, "main", "png")

RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")

RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
DAT <- "04_Figures/F09/c_data"
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED  <- unname(DIR_COLORS["Up"])
COMP_BLUE <- unname(DIR_COLORS["Down"])
N_SHOW <- 5

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                           show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")

scatter_df <- dep_df %>%
  transmute(gene, logFC_HR = logFC_Acute_HR, logFC_LR = logFC_Acute_LR,
            pi_HR = pi_score_Acute_HR, pi_LR = pi_score_Acute_LR) %>%
  filter(!is.na(logFC_HR), !is.na(logFC_LR)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(imputed = replace_na(imputed, FALSE),
         sig_class = classify_concordance(pi_HR, pi_LR, "HR", "LR"),
         sig_class = factor(sig_class, levels = names(SIG_COLORS_F07)),
         is_sig = sig_class != "NS",
         quadrant = case_when(
           logFC_HR > 0 & logFC_LR > 0 ~ "Concordant Up",
           logFC_HR < 0 & logFC_LR < 0 ~ "Concordant Down",
           logFC_HR > 0 & logFC_LR < 0 ~ "Discordant HR Up",
           TRUE ~ "Discordant LR Up"))

universe <- scatter_df$gene
message(sprintf("  Total: %d | Sig: %d", nrow(scatter_df), sum(scatter_df$is_sig)))

pw_collection <- build_pathway_collection(min_size = 15, max_size = 500)

run_set_ora <- function(genes, set_name) {
  if (length(genes) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(genes = genes, universe = universe, pathways = pw_collection,
                          jaccard_cutoff = 0.5, min_size = 15, max_size = 500, padj_cutoff = 1),
    error = function(e) { message("  ORA error: ", e$message); tibble() })
  if (nrow(res) == 0) return(tibble())
  res %>% mutate(set = set_name, pathway_label = clean_pathway_name(pathway),
                  neg_log10_padj = -log10(padj), significant = padj < 0.05) %>%
    arrange(desc(neg_log10_padj)) %>% slice_head(n = N_SHOW)
}

ora_q1 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Concordant Up"], "Concordant Up")
ora_q2 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Discordant LR Up"], "Discordant LR Up")
ora_q3 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Concordant Down"], "Concordant Down")
ora_q4 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Discordant HR Up"], "Discordant HR Up")
all_quad_ora <- bind_rows(ora_q1, ora_q2, ora_q3, ora_q4)
if (nrow(all_quad_ora) > 0)
  write_csv(all_quad_ora, file.path(DAT, "panel_A", "ora_quadrant.csv"))

# Scatter
xlim_range <- range(scatter_df$logFC_HR, na.rm = TRUE) * 1.15
ylim_range <- range(scatter_df$logFC_LR, na.rm = TRUE) * 1.15
ns_df  <- filter(scatter_df, sig_class == "NS")
sig_df <- filter(scatter_df, sig_class != "NS")

q_df <- scatter_df %>% mutate(q = case_when(
  logFC_HR > 0 & logFC_LR > 0 ~ "Q1", logFC_HR < 0 & logFC_LR < 0 ~ "Q3",
  logFC_HR > 0 & logFC_LR < 0 ~ "Q4", TRUE ~ "Q2"))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig <- q_df %>% filter(is_sig) %>% count(q) %>% deframe()
for (qq in c("Q1","Q2","Q3","Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- sig_df %>% group_by(sig_class) %>%
  arrange(desc(abs(logFC_HR) + abs(logFC_LR))) %>% slice_head(n = 5) %>% ungroup() %>%
  mutate(label_fill = SIG_LABEL_FILL_F07[as.character(sig_class)],
         label_text_col = SIG_LABEL_TEXT_F07[as.character(sig_class)])

txt_gene <- scale_text(BASE_GENE, 250) * 1.2
txt_quad <- scale_text(BASE_QUADRANT, 250) * 1.3

p_scatter <- ggplot(mapping = aes(x = logFC_HR, y = logFC_LR)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, color = "grey80", fill = "grey85", shape = 21,
             size = 1.0, alpha = 0.3, stroke = 0.2) +
  geom_point(data = sig_df, aes(fill = sig_class), shape = 21, size = 2.3,
             color = ifelse(sig_df$imputed, "black", "grey75"), alpha = 0.80, stroke = 0.9) +
  scale_fill_manual(values = SIG_COLORS_F07, name = "Significance") +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   size = txt_gene, fontface = "italic", max.overlaps = 50,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.4, point.padding = 0.3, force = 4, force_pull = 0.3,
                   label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  annotate("label", x = xlim_range[2], y = ylim_range[2],
           label = sprintf("Concordant Up  %s/%s", q_sig["Q1"], q_counts["Q1"]),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92), label.padding = unit(2.5, "pt")) +
  annotate("label", x = xlim_range[1], y = ylim_range[1],
           label = sprintf("Concordant Down  %s/%s", q_sig["Q3"], q_counts["Q3"]),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92), label.padding = unit(2.5, "pt")) +
  annotate("label", x = xlim_range[1], y = ylim_range[2],
           label = sprintf("Discordant (LR Up)  %s/%s", q_sig["Q2"], q_counts["Q2"]),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92), label.padding = unit(2.5, "pt")) +
  annotate("label", x = xlim_range[2], y = ylim_range[1],
           label = sprintf("Discordant (HR Up)  %s/%s", q_sig["Q4"], q_counts["Q4"]),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92), label.padding = unit(2.5, "pt")) +
  annotate("text", x = mean(xlim_range), y = ylim_range[1],
           label = "log2FC (Acute HR)", hjust = 0.5, vjust = -0.4,
           size = 3.2, color = "grey30", fontface = "bold") +
  annotate("text", x = xlim_range[1], y = mean(ylim_range),
           label = "log2FC (Acute LR)", hjust = 0.5, vjust = -0.4,
           size = 3.2, color = "grey30", fontface = "bold", angle = 90) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(x = NULL, y = NULL) + FIG_THEME +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), plot.margin = margin(2, 0, 2, 0, "mm"),
        legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10), legend.key.size = unit(5, "mm"),
        legend.margin = margin(-5, 0, 0, 0), legend.box.margin = margin(-8, 0, 0, 0)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 5, alpha = 0.8)))

n_total <- nrow(scatter_df); n_sig <- sum(scatter_df$is_sig)
r_pear <- cor(scatter_df$logFC_HR, scatter_df$logFC_LR, use = "complete.obs")

p_final <- p_scatter + plot_annotation(
  title = "Acute Concordance: Quadrant ORA",
  subtitle = sprintf("N = %d | %d DEPs (Pi < 0.05) | r = %.2f | Acute_HR (x) vs Acute_LR (y)",
                      n_total, n_sig, r_pear),
  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey30")))

ggsave(file.path(RPT_SUPP_PDF, "c01_ORA_SUPP.pdf"), p_final,
       width = 250, height = 250, units = "mm", device = pdf_device)
ggsave(file.path(RPT_SUPP_PNG, "c01_ORA_SUPP.png"), p_final,
       width = 250, height = 250, units = "mm", dpi = 300)
message("F09 Supp C01 ORA done")
