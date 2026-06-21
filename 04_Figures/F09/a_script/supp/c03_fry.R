# F09 Supp C03: fry Rotation Test — Acute Concordance
# Acute_HR-sig sets tested in Acute_LR. (Moved from main to supplementary)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F09/a_script/style.R")
suppressPackageStartupMessages({
  library(tidyverse); library(limma); library(fgsea); library(patchwork)
})
set.seed(42)

RPT <- "04_Figures/F09/b_reports"
RPT_PDF       <- file.path(RPT, "main", "pdf")
RPT_PNG       <- file.path(RPT, "main", "png")
RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")
RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
DAT <- "04_Figures/F09/c_data"
dir.create(file.path(DAT, "panel_C_fry"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device(); PF_W <- 220

dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
imp_csv <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)

meta <- dal$metadata; sample_cols <- meta$Col_ID
mat_imp <- imp_csv %>% select(uniprot_id, all_of(sample_cols)) %>%
  column_to_rownames("uniprot_id") %>% as.matrix()

meta$Group_Time <- factor(meta$Group_Time,
  levels = c("HR_T1","HR_T2","HR_T3","LR_T1","LR_T2","LR_T3"))
design <- model.matrix(~ 0 + Group_Time, data = meta)
colnames(design) <- gsub("^Group_Time", "", colnames(design))
block_id <- sub("_T[123]$", "", meta$Col_ID)
corfit <- duplicateCorrelation(mat_imp, design, block = block_id)
cor_imp <- corfit$consensus.correlation

cm <- makeContrasts(Acute_LR = LR_T3 - LR_T2, levels = design)

imp_ids <- rownames(mat_imp)
sig_pi <- dep_df %>% filter(pi_score_Acute_HR < 0.05, uniprot_id %in% imp_ids)
sets_pi <- list(
  up = match(sig_pi$uniprot_id[sig_pi$logFC_Acute_HR > 0], imp_ids),
  down = match(sig_pi$uniprot_id[sig_pi$logFC_Acute_HR < 0], imp_ids),
  up_ids = sig_pi$uniprot_id[sig_pi$logFC_Acute_HR > 0],
  down_ids = sig_pi$uniprot_id[sig_pi$logFC_Acute_HR < 0])
message(sprintf("Gene sets: up=%d, down=%d", length(sets_pi$up), length(sets_pi$down)))

run_fry_set <- function(idx, nm) {
  if (length(idx) < 3) return(tibble(set=nm, n=length(idx), direction=NA, PValue=NA_real_))
  res <- fry(mat_imp, index=idx, design=design, contrast=cm[,"Acute_LR"],
             block=block_id, correlation=cor_imp)
  tibble(set=nm, n=length(idx), direction=res$Direction[1],
         PValue=res$PValue[1], PValue.Mixed=res$PValue.Mixed[1])
}

fry_up <- run_fry_set(sets_pi$up, "ac_hr_up") %>% mutate(expected="Up", consistent=direction==expected)
fry_dn <- run_fry_set(sets_pi$down, "ac_hr_down") %>% mutate(expected="Down", consistent=direction==expected)
fry_all <- bind_rows(fry_up, fry_dn) %>% mutate(cor_within=cor_imp)
write_csv(fry_all, file.path(DAT, "panel_C_fry", "fry_results.csv"))

t_rank <- dep_df %>% filter(uniprot_id %in% imp_ids, !is.na(t_Acute_LR)) %>%
  arrange(desc(t_Acute_LR)) %>%
  mutate(rank=row_number(), in_up=uniprot_id %in% sets_pi$up_ids,
         in_down=uniprot_id %in% sets_pi$down_ids)

running_es <- function(t_vals, in_set) {
  n <- length(t_vals); n_h <- sum(in_set)
  if (n_h == 0) return(rep(0, n))
  hit_w <- ifelse(in_set, abs(t_vals), 0); miss_w <- 1/(n-n_h)
  cumsum(ifelse(in_set, hit_w/sum(hit_w), -miss_w))
}
t_rank$es_up <- running_es(t_rank$t_Acute_LR, t_rank$in_up)
t_rank$es_down <- running_es(t_rank$t_Acute_LR, t_rank$in_down)
n_all <- nrow(t_rank); txt_s <- scale_text(BASE_STAT, PF_W)

HR_COL <- unname(CONTRAST_COLORS["Acute_HR"])
LR_COL <- unname(CONTRAST_COLORS["Acute_LR"])

make_bc <- function(t_df, in_col, es_col, fr, title, col) {
  marks <- t_df %>% filter(.data[[in_col]])
  is_sig <- !is.na(fr$PValue) && fr$PValue < 0.05
  lc <- if (is_sig) col else "grey55"
  p_es <- ggplot(t_df, aes(x=rank, y=.data[[es_col]])) +
    geom_area(fill=scales::alpha(lc, 0.15)) + geom_line(color=lc, linewidth=0.6) +
    geom_hline(yintercept=0, linetype="dashed", color="grey60", linewidth=0.3) +
    annotate("text", x=n_all*0.98, y=Inf,
             label=sprintf("fry %s, %s (n=%d)", fr$direction, fmt_p(fr$PValue), fr$n),
             hjust=1, vjust=1.3, size=txt_s*1.1, fontface="bold",
             color=if(fr$consistent) "grey20" else "#DC2626") +
    labs(y="ES", title=title) +
    scale_x_continuous(limits=c(1,n_all), expand=c(0.005,0)) + FIG_THEME +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), plot.margin=margin(3,4,0,4,"mm"),
          plot.title=element_text(size=9.5, face="bold"))
  p_bc <- ggplot(marks, aes(x=rank, xend=rank, y=0, yend=1)) +
    geom_segment(color=lc, linewidth=0.3, alpha=0.7) +
    scale_x_continuous(limits=c(1,n_all), expand=c(0.005,0)) +
    scale_y_continuous(expand=c(0,0)) + FIG_THEME +
    theme(axis.text=element_blank(), axis.title=element_blank(),
          axis.ticks=element_blank(), panel.grid=element_blank(),
          panel.background=element_rect(fill="grey97"), plot.margin=margin(0,4,0,4,"mm"))
  list(es=p_es, bc=p_bc)
}

p1 <- make_bc(t_rank, "in_up", "es_up", fry_up,
              sprintf("Ac.(HR)-Up DEPs (n=%d) -> Ac.(LR) ranked t", length(sets_pi$up)), HR_COL)
p2 <- make_bc(t_rank, "in_down", "es_down", fry_dn,
              sprintf("Ac.(HR)-Down DEPs (n=%d) -> Ac.(LR) ranked t%s",
                      length(sets_pi$down),
                      if (!is.na(fry_dn$PValue) && fry_dn$PValue > 0.05) "  (n.s.)" else ""), HR_COL)

p_t <- ggplot(t_rank, aes(x=rank, y=t_Acute_LR)) +
  geom_area(fill=scales::alpha(LR_COL, 0.20), color=LR_COL, linewidth=0.3) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  labs(x=sprintf("Protein rank by t(Acute LR) [n=%d]", n_all), y="t-stat") +
  scale_x_continuous(limits=c(1,n_all), expand=c(0.005,0)) +
  FIG_THEME + theme(plot.margin=margin(2,4,4,4,"mm"))

pC <- p1$es / p1$bc / p2$es / p2$bc / p_t +
  plot_layout(heights=c(2.5,0.4,2.5,0.4,1.2)) +
  plot_annotation(title="fry: Acute Concordance (HR-sig in LR)",
    subtitle=sprintf("dupCor=%.3f | n=%d | Acute_HR-sig sets tested in Acute_LR", cor_imp, n_all),
    theme=theme(plot.title=element_text(size=11, face="bold"),
                plot.subtitle=element_text(size=8.5, color="grey30")))

ggsave(file.path(RPT_SUPP_PDF, "c03_fry_SUPP.pdf"), pC,
       width=PF_W, height=175, units="mm", device=pdf_device)
ggsave(file.path(RPT_SUPP_PNG, "c03_fry_SUPP.png"), pC,
       width=PF_W, height=175, units="mm", dpi=300)
message("F09 Supp C03 fry done")
