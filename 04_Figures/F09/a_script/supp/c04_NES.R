# F09 Supp C04: fGSEA NES Scatter — Acute Concordance (Acute_HR vs Acute_LR)
# (Moved from main panel to supplementary)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F09/a_script/style.R")
suppressPackageStartupMessages({ library(tidyverse); library(fgsea); library(ggrepel) })

PG_W <- 200
RPT <- "04_Figures/F09/b_reports"
RPT_PDF       <- file.path(RPT, "main", "pdf")
RPT_PNG       <- file.path(RPT, "main", "png")
RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")
RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
DAT <- "04_Figures/F09/c_data"
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/a_non_imputed/c_data/03_combined_results.csv", show_col_types = FALSE)
pw_collection <- build_pathway_collection(min_size = 10, max_size = 500)

set.seed(42)
fgsea_list <- list()
for (ctr in c("Acute_HR", "Acute_LR")) {
  stats <- dep_df[[paste0("t_", ctr)]]
  names(stats) <- dep_df$gene
  stats <- sort(stats[!is.na(stats) & is.finite(stats)], decreasing = TRUE)
  res <- run_fgsea_deduplicated(ranks = stats, pathways = pw_collection,
                                 jaccard_cutoff = 0.5, nperm = 10000,
                                 min_size = 10, max_size = 500)
  res$contrast <- ctr
  fgsea_list[[ctr]] <- res
  message(sprintf("  fgsea done: %s (%d terms)", ctr, nrow(res)))
}
fgsea_all <- bind_rows(fgsea_list)

avail_dbs <- unique(fgsea_all$database)
use_dbs <- intersect(c("Hallmark", "GO:BP"), avail_dbs)
if (length(use_dbs) == 0) use_dbs <- avail_dbs[1:min(2, length(avail_dbs))]

fgsea_wide <- fgsea_all %>%
  filter(database %in% use_dbs, contrast %in% c("Acute_HR", "Acute_LR")) %>%
  select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Acute_HR), !is.na(NES_Acute_LR)) %>%
  mutate(set_size = coalesce(size_Acute_HR, size_Acute_LR),
         sig_HR = !is.na(padj_Acute_HR) & padj_Acute_HR < 0.05,
         sig_LR = !is.na(padj_Acute_LR) & padj_Acute_LR < 0.05,
         significance = case_when(
           sig_HR & sig_LR ~ "Sig Both",
           sig_HR          ~ "Sig HR only",
           sig_LR          ~ "Sig LR only",
           TRUE            ~ "NS"
         ) %>% factor(levels = names(SIG_COLORS_F07)),
         pathway_label = clean_pathway_name(pathway))

fgsea_sig <- filter(fgsea_wide, significance != "NS")
message(sprintf("  %d pathways | %d sig", nrow(fgsea_wide), nrow(fgsea_sig)))

nes_cor <- cor.test(fgsea_wide$NES_Acute_HR, fgsea_wide$NES_Acute_LR, method = "spearman")
nes_ci <- fisher_z_ci(nes_cor$estimate, nrow(fgsea_wide))
nes_lim <- max(abs(c(fgsea_wide$NES_Acute_HR, fgsea_wide$NES_Acute_LR))) * 1.15

n_conc_up <- sum(fgsea_sig$NES_Acute_HR > 0 & fgsea_sig$NES_Acute_LR > 0)
n_conc_dn <- sum(fgsea_sig$NES_Acute_HR < 0 & fgsea_sig$NES_Acute_LR < 0)

label_pw <- fgsea_sig %>%
  mutate(label_fill = SIG_LABEL_FILL_F07[as.character(significance)],
         label_text_col = SIG_LABEL_TEXT_F07[as.character(significance)])
ns_df <- filter(fgsea_wide, significance == "NS")
sig_df <- filter(fgsea_wide, significance != "NS")

txt_pw <- scale_text(BASE_PATHWAY, PG_W) * 1.25
txt_quad <- scale_text(BASE_QUADRANT, PG_W) * 1.4

pD <- ggplot(mapping = aes(x = NES_Acute_HR, y = NES_Acute_LR)) +
  annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf,
           fill="#FFE0E0", alpha=0.55, color="grey70", linewidth=0.2) +
  annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=0,
           fill="#FFE0E0", alpha=0.55, color="grey70", linewidth=0.2) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0,
           fill="#DCEEFF", alpha=0.55, color="grey70", linewidth=0.2) +
  annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf,
           fill="#DCEEFF", alpha=0.55, color="grey70", linewidth=0.2) +
  geom_hline(yintercept=0, color="grey60", linewidth=0.2) +
  geom_vline(xintercept=0, color="grey60", linewidth=0.2) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="black", linewidth=0.3) +
  geom_point(data=ns_df, aes(shape=database),
             size=2.5, fill="grey70", color="grey55", alpha=0.40, stroke=0.4) +
  geom_point(data=sig_df, aes(fill=significance, size=set_size, shape=database),
             color="grey65", alpha=0.80, stroke=0.8) +
  scale_fill_manual(values=SIG_COLORS_F07, name="Significance") +
  scale_shape_manual(values=c("Hallmark"=24, "GO:BP"=21), name="Database") +
  scale_size_continuous(range=c(3,10), name="Set size", breaks=c(20,50,100,200)) +
  geom_label_repel(data=label_pw, aes(label=pathway_label),
                   fill=label_pw$label_fill, color=label_pw$label_text_col,
                   size=txt_pw, fontface="bold", max.overlaps=50,
                   segment.size=0.2, segment.color="grey50",
                   min.segment.length=0, show.legend=FALSE,
                   box.padding=0.5, point.padding=0.4, force=5, force_pull=0.3,
                   label.padding=unit(1.5,"pt"), label.r=unit(1,"pt"),
                   label.size=0.15, seed=42) +
  annotate("label", x=Inf, y=Inf, label=sprintf("Concordant Up  n = %d", n_conc_up),
           hjust=1.05, vjust=1.2, size=txt_quad, fontface="bold",
           color="#D6604D", fill=alpha("white",0.92), label.padding=unit(2.5,"pt")) +
  annotate("label", x=-Inf, y=-Inf, label=sprintf("Concordant Down  n = %d", n_conc_dn),
           hjust=-0.05, vjust=-0.2, size=txt_quad, fontface="bold",
           color="#D6604D", fill=alpha("white",0.92), label.padding=unit(2.5,"pt")) +
  coord_cartesian(xlim=c(-nes_lim, nes_lim), ylim=c(-nes_lim, nes_lim)) +
  labs(title="Pathway Concordance: Acute (HR vs LR)",
       subtitle=sprintf("%s | %d pathways (%d sig.) | rho = %.2f [%.2f, %.2f]",
                         paste(use_dbs, collapse=" + "), nrow(fgsea_wide), nrow(fgsea_sig),
                         nes_cor$estimate, nes_ci[1], nes_ci[2]),
       x="NES (Acute HR)", y="NES (Acute LR)") +
  FIG_THEME +
  theme(legend.position="bottom", legend.box="horizontal") +
  guides(fill=guide_legend(nrow=1, order=1, override.aes=list(size=4, alpha=0.8, shape=21)),
         shape=guide_legend(nrow=1, order=2, override.aes=list(size=4, fill="grey50")),
         size=guide_legend(nrow=1, order=3))

ggsave(file.path(RPT_SUPP_PDF, "c04_NES_SUPP.pdf"), pD,
       width=PG_W, height=PG_W, units="mm", device=pdf_device)
ggsave(file.path(RPT_SUPP_PNG, "c04_NES_SUPP.png"), pD,
       width=PG_W, height=PG_W, units="mm", dpi=300)
message("F09 Supp C04 NES done")
