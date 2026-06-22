# F06 Panel C (LR): fry Rotation Test — Training-Sig Tested in Acute
# Defines gene sets from Training_LR Pi<0.05, tests against Acute_LR contrast.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(fgsea)
  library(patchwork)
})
set.seed(42)

GRP <- "LR"
CTR_SIG <- "Training_LR"
CTR_TEST <- "Acute_LR"
RPT <- "04_Figures/F06/b_reports"
DAT <- "04_Figures/F06/c_data"
dir.create(file.path(DAT, paste0("panel_C_", GRP)), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()
PF_W <- 220

dal <- readRDS("02_Imputation/c_data/DAList_imputed_missforest.rds")
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
imp_csv <- dplyr::bind_cols(
  tibble::as_tibble(dal$annotation[, c("uniprot_id", "protein", "gene", "description")]),
  tibble::as_tibble(dal$data)
)

meta <- dal$metadata
sample_cols <- meta$Col_ID

mat_imp <- imp_csv %>%
  select(uniprot_id, all_of(sample_cols)) %>%
  column_to_rownames("uniprot_id") %>%
  as.matrix()

meta$Group_Time <- factor(meta$Group_Time,
  levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
)
design <- model.matrix(~ 0 + Group_Time, data = meta)
colnames(design) <- gsub("^Group_Time", "", colnames(design))

block_id <- sub("_T[123]$", "", meta$Col_ID)
corfit_imp <- duplicateCorrelation(mat_imp, design, block = block_id)
cor_imp <- corfit_imp$consensus.correlation

cm <- makeContrasts(Acute_LR = LR_T3 - LR_T2, levels = design)

imp_ids <- rownames(mat_imp)
sig_pi <- dep_df %>%
  filter(.data[[paste0("pi_score_", CTR_SIG)]] < 0.05, uniprot_id %in% imp_ids)

sets_pi <- list(
  up       = match(sig_pi$uniprot_id[sig_pi[[paste0("logFC_", CTR_SIG)]] > 0], imp_ids),
  down     = match(sig_pi$uniprot_id[sig_pi[[paste0("logFC_", CTR_SIG)]] < 0], imp_ids),
  up_ids   = sig_pi$uniprot_id[sig_pi[[paste0("logFC_", CTR_SIG)]] > 0],
  down_ids = sig_pi$uniprot_id[sig_pi[[paste0("logFC_", CTR_SIG)]] < 0]
)
message(sprintf("Gene sets (Pi < 0.05): up = %d, down = %d", length(sets_pi$up), length(sets_pi$down)))

run_fry_set <- function(idx, set_name) {
  if (length(idx) < 3) {
    return(tibble(
      set = set_name, n = length(idx),
      direction = NA_character_, PValue = NA_real_
    ))
  }
  res <- fry(mat_imp,
    index = idx, design = design,
    contrast = cm[, "Acute_LR"], block = block_id, correlation = cor_imp
  )
  tibble(
    set = set_name, n = length(idx), direction = res$Direction[1],
    PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1]
  )
}

fry_up <- run_fry_set(sets_pi$up, paste0(GRP, "_up")) %>%
  mutate(expected = "Up", consistent = direction == expected)
fry_dn <- run_fry_set(sets_pi$down, paste0(GRP, "_down")) %>%
  mutate(expected = "Down", consistent = direction == expected)
fry_all <- bind_rows(fry_up, fry_dn) %>% mutate(cor_within = cor_imp)
write_csv(fry_all, file.path(DAT, paste0("panel_C_", GRP), "fry_results.csv"))

# Barcode
t_col <- paste0("t_", CTR_TEST)
t_rank <- dep_df %>%
  filter(uniprot_id %in% imp_ids, !is.na(.data[[t_col]])) %>%
  arrange(desc(.data[[t_col]])) %>%
  mutate(
    rank = row_number(),
    in_up = uniprot_id %in% sets_pi$up_ids,
    in_down = uniprot_id %in% sets_pi$down_ids
  )

running_es <- function(t_vals, in_set) {
  n <- length(t_vals)
  n_h <- sum(in_set)
  if (n_h == 0) {
    return(rep(0, n))
  }
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  miss_w <- 1 / (n - n_h)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -miss_w))
}

t_rank$es_up <- running_es(t_rank[[t_col]], t_rank$in_up)
t_rank$es_down <- running_es(t_rank[[t_col]], t_rank$in_down)
n_all <- nrow(t_rank)
txt_s <- scale_text(BASE_STAT, PF_W)

TR_COLOR <- unname(CONTRAST_COLORS[CTR_SIG])
AC_COLOR <- unname(CONTRAST_COLORS[CTR_TEST])

make_barcode <- function(t_df, in_col, es_col, fry_row, title, color) {
  marks <- t_df %>% filter(.data[[in_col]])
  is_sig <- !is.na(fry_row$PValue) && fry_row$PValue < 0.05
  line_color <- if (is_sig) color else "grey55"
  p_es <- ggplot(t_df, aes(x = rank, y = .data[[es_col]])) +
    geom_area(fill = scales::alpha(line_color, 0.15)) +
    geom_line(color = line_color, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    annotate("text",
      x = n_all * 0.98, y = Inf,
      label = sprintf(
        "fry %s, %s (n = %d)", fry_row$direction,
        fmt_p(fry_row$PValue), fry_row$n
      ),
      hjust = 1, vjust = 1.3, size = txt_s * 1.1, fontface = "bold",
      color = if (fry_row$consistent) "grey20" else "#DC2626"
    ) +
    labs(y = "ES", title = title) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      axis.ticks.x = element_blank(), plot.margin = margin(3, 4, 0, 4, "mm"),
      plot.title = element_text(size = 9.5, face = "bold")
    )
  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(
      axis.text = element_blank(), axis.title = element_blank(),
      axis.ticks = element_blank(), panel.grid = element_blank(),
      panel.background = element_rect(fill = "grey97"),
      plot.margin = margin(0, 4, 0, 4, "mm")
    )
  list(es = p_es, bc = p_bc)
}

up_title <- sprintf(
  "Tr.(%s)-Up DEPs (n = %d) -> %s ranked t",
  GRP, length(sets_pi$up), CTR_SHORT[CTR_TEST]
)
dn_title <- sprintf(
  "Tr.(%s)-Down DEPs (n = %d) -> %s ranked t%s",
  GRP, length(sets_pi$down), CTR_SHORT[CTR_TEST],
  if (!is.na(fry_dn$PValue) && fry_dn$PValue > 0.05) "  (n.s.)" else ""
)

p1 <- make_barcode(t_rank, "in_up", "es_up", fry_up, up_title, TR_COLOR)
p2 <- make_barcode(t_rank, "in_down", "es_down", fry_dn, dn_title, TR_COLOR)

p_t <- ggplot(t_rank, aes(x = rank, y = .data[[t_col]])) +
  geom_area(fill = scales::alpha(AC_COLOR, 0.20), color = AC_COLOR, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = sprintf("Protein rank by t(%s)  [n = %d]", CTR_FACET[CTR_TEST], n_all), y = "t-stat") +
  scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(plot.margin = margin(2, 4, 4, 4, "mm"))

pC <- p1$es / p1$bc / p2$es / p2$bc / p_t +
  plot_layout(heights = c(2.5, 0.4, 2.5, 0.4, 1.2)) +
  plot_annotation(
    title = sprintf("fry Rotation Test: Training-Sig in Acute (%s)", GRP),
    subtitle = sprintf(
      "dupCor = %.3f | n = %d proteins | %s-sig sets tested in %s",
      cor_imp, n_all, CTR_SHORT[CTR_SIG], CTR_SHORT[CTR_TEST]
    ),
    theme = theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "grey30")
    )
  )

ggsave(file.path(RPT, paste0("panel_C_fry_", GRP, "_MAIN.pdf")), pC,
  width = PF_W, height = 175, units = "mm", device = pdf_device
)
ggsave(file.path(RPT, paste0("panel_C_fry_", GRP, "_MAIN.png")), pC,
  width = PF_W, height = 175, units = "mm", dpi = 300
)

message(sprintf("F06 Panel C fry (%s) done", GRP))
