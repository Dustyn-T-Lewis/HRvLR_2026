# F05 Panel D: fry Rotation Test — Between-Group Concordance Barcode
# Tests whether HR Pi-sig proteins respond concordantly in LR contrast,
# and vice versa, for both Training and Acute pairs.
#
# Training pair: Training_HR Pi-sig sets tested against Training_LR t-stat ranks
# Acute pair:   Acute_HR Pi-sig sets tested against Acute_LR t-stat ranks
#
# Design: HR and LR are independent subject groups -> no circularity.
# Reference: Wu & Smyth 2010, Bioinformatics -- ROAST/fry
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F05/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(patchwork)
})

set.seed(42)

PD_W <- 400
PD_H <- 150

RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(file.path(DAT, "panel_D_fry"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# Load data

dal <- readRDS("02_Imputation/c_data/DAList_imputed_missforest.rds")
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
meta <- read_csv("00_input/HRvLR_meta.csv", show_col_types = FALSE)

imp_csv <- dplyr::bind_cols(
  tibble::as_tibble(dal$annotation[, c("uniprot_id", "protein", "gene", "description")]),
  tibble::as_tibble(dal$data)
)

# Intersect metadata with actual imputed samples (outliers removed during QC)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
avail_samples <- setdiff(names(imp_csv), ann_cols)
sample_cols <- intersect(meta$Col_ID, avail_samples)
meta <- meta %>% filter(Col_ID %in% sample_cols)

# Imputed matrix

mat_imp <- imp_csv %>%
  select(uniprot_id, all_of(sample_cols)) %>%
  column_to_rownames("uniprot_id") %>%
  as.matrix()

n_imp <- nrow(mat_imp)
message(sprintf("Imputed matrix: %d proteins x %d samples", n_imp, ncol(mat_imp)))

# Design + duplicateCorrelation

meta$subject <- sub("_T[123]$", "", meta$Col_ID)
meta$Group_Time <- factor(meta$Group_Time)
design <- model.matrix(~ 0 + Group_Time, data = meta)
colnames(design) <- gsub("^Group_Time", "", colnames(design))

block_id <- meta$subject
corfit <- duplicateCorrelation(mat_imp, design, block = block_id)
cor_val <- corfit$consensus.correlation
message(sprintf("Within-subject correlation: %.4f", cor_val))

# Available design columns
avail <- colnames(design)
message(sprintf("Design columns: %s", paste(avail, collapse = ", ")))

# fry-based barcode for both pairs

# Pair definitions: source contrast (define gene sets) -> target contrast (rank against)
pair_defs <- list(
  Training = list(
    src_lfc = "logFC_Training_HR", src_pi = "pi_score_Training_HR",
    tgt_t = "t_Training_LR",
    src_label = "Training HR", tgt_label = "Training LR",
    # Target contrast: LR_T2 - LR_T1
    tgt_t2_pat = "LR.*T2", tgt_t1_pat = "LR.*T1"
  ),
  Acute = list(
    src_lfc = "logFC_Acute_HR", src_pi = "pi_score_Acute_HR",
    tgt_t = "t_Acute_LR",
    src_label = "Acute HR", tgt_label = "Acute LR",
    # Target contrast: LR_T3 - LR_T2
    tgt_t2_pat = "LR.*T3", tgt_t1_pat = "LR.*T2"
  )
)

imp_ids <- rownames(mat_imp)

# Running enrichment score
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

fry_results_all <- list()
barcode_plots <- list()

for (pair_name in names(pair_defs)) {
  pd <- pair_defs[[pair_name]]

  # Build contrast matrix for target
  t2_cols <- avail[grepl(pd$tgt_t2_pat, avail)]
  t1_cols <- avail[grepl(pd$tgt_t1_pat, avail)]

  if (length(t2_cols) == 0 || length(t1_cols) == 0) {
    message(sprintf("  Skipping %s fry -- cannot build target contrast", pair_name))
    next
  }

  cm <- matrix(0, nrow = ncol(design), ncol = 1)
  rownames(cm) <- colnames(design)
  colnames(cm) <- "target"
  for (col in t2_cols) cm[col, 1] <- 1 / length(t2_cols)
  for (col in t1_cols) cm[col, 1] <- -1 / length(t1_cols)

  # Define source gene sets (Pi < 0.05)
  sig_df <- dep_df %>%
    filter(
      !is.na(.data[[pd$src_pi]]) & .data[[pd$src_pi]] < 0.05,
      uniprot_id %in% imp_ids
    )

  idx_up <- match(sig_df$uniprot_id[sig_df[[pd$src_lfc]] > 0], imp_ids)
  idx_down <- match(sig_df$uniprot_id[sig_df[[pd$src_lfc]] < 0], imp_ids)
  idx_up <- idx_up[!is.na(idx_up)]
  idx_down <- idx_down[!is.na(idx_down)]

  message(sprintf(
    "  %s: %d up, %d down Pi-sig proteins", pair_name,
    length(idx_up), length(idx_down)
  ))

  # Run fry
  fry_list <- list()
  for (dir in c("up", "down")) {
    idx <- if (dir == "up") idx_up else idx_down
    if (length(idx) < 3) {
      fry_list[[dir]] <- tibble(
        pair = pair_name, set = dir, n = length(idx),
        direction = NA_character_, PValue = NA_real_,
        PValue.Mixed = NA_real_
      )
      next
    }
    res <- fry(mat_imp,
      index = idx, design = design,
      contrast = cm[, "target"], block = block_id, correlation = cor_val
    )
    fry_list[[dir]] <- tibble(
      pair = pair_name, set = dir, n = length(idx),
      direction = res$Direction[1],
      PValue = res$PValue[1],
      PValue.Mixed = res$PValue.Mixed[1]
    )
  }

  fry_pair <- bind_rows(fry_list) %>%
    mutate(
      expected = ifelse(set == "up", "Up", "Down"),
      consistent = direction == expected
    )

  fry_results_all[[pair_name]] <- fry_pair

  # Barcode data
  up_ids <- sig_df$uniprot_id[sig_df[[pd$src_lfc]] > 0]
  down_ids <- sig_df$uniprot_id[sig_df[[pd$src_lfc]] < 0]

  t_rank <- dep_df %>%
    filter(uniprot_id %in% imp_ids, !is.na(.data[[pd$tgt_t]])) %>%
    arrange(desc(.data[[pd$tgt_t]])) %>%
    mutate(
      rank = row_number(),
      in_up = uniprot_id %in% up_ids,
      in_down = uniprot_id %in% down_ids
    )

  t_rank$es_up <- running_es(t_rank[[pd$tgt_t]], t_rank$in_up)
  t_rank$es_down <- running_es(t_rank[[pd$tgt_t]], t_rank$in_down)

  fry_up <- fry_pair %>% filter(set == "up")
  fry_dn <- fry_pair %>% filter(set == "down")
  n_all <- nrow(t_rank)
  txt_s <- scale_text(BASE_STAT, PD_W / 2)

  # Up barcode + ES
  marks_up <- t_rank %>% filter(in_up)

  p_es_up <- ggplot(t_rank, aes(x = rank, y = es_up)) +
    annotate("rect",
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = "#D6604D", alpha = 0.04
    ) +
    geom_area(fill = alpha("#D6604D", 0.2), color = NA) +
    geom_line(color = "#D6604D", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    annotate("text",
      x = n_all * 0.98, y = Inf,
      label = sprintf(
        "fry %s, %s (n=%d)",
        fry_up$direction, fmt_p(fry_up$PValue), fry_up$n
      ),
      hjust = 1, vjust = 1.3, size = txt_s * 1.1, fontface = "bold",
      color = ifelse(isTRUE(fry_up$consistent), "grey20", "#DC2626")
    ) +
    labs(y = "ES", title = sprintf(
      "%s-Up DEPs -> %s t-stats",
      pd$src_label, pd$tgt_label
    )) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(4, 4, 0, 4, "mm"),
      plot.title = element_text(size = 9, face = "bold")
    )

  p_bc_up <- ggplot(marks_up, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = "#D6604D", linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(
      axis.text = element_blank(), axis.title = element_blank(),
      axis.ticks = element_blank(), panel.grid = element_blank(),
      panel.background = element_rect(fill = "grey97"),
      plot.margin = margin(0, 4, 0, 4, "mm")
    )

  # t-stat area
  p_t <- ggplot(t_rank, aes(x = rank, y = .data[[pd$tgt_t]])) +
    geom_area(fill = alpha("#5DA5DA", 0.25), color = "#5DA5DA", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    labs(
      x = sprintf("Rank by t(%s)  [%d proteins]", pd$tgt_label, n_all),
      y = "t-stat"
    ) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(plot.margin = margin(2, 4, 4, 4, "mm"))

  barcode_plots[[pair_name]] <- p_es_up / p_bc_up / p_t +
    plot_layout(heights = c(3, 0.5, 1.5))
}

fry_all <- bind_rows(fry_results_all)
write_csv(fry_all, file.path(DAT, "panel_D_fry", "fry_results.csv"))

# -- Composite: side-by-side barcode panels ------------------------------------

if (length(barcode_plots) == 2) {
  pD <- (barcode_plots$Training | barcode_plots$Acute) +
    plot_annotation(
      title = "fry Rotation Test: HR DEP sets against LR ranked t-statistics",
      subtitle = sprintf(
        "missForest-imputed | dupCor=%.3f | Pi<0.05 gene sets",
        cor_val
      ),
      tag_levels = list(c("D", "")),
      theme = theme(
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, face = "bold.italic", color = "grey30"),
        plot.tag      = element_text(size = 15, face = "bold")
      )
    )
} else if (length(barcode_plots) == 1) {
  pD <- barcode_plots[[1]] +
    plot_annotation(
      tag_levels = list(c("D")),
      theme = theme(plot.tag = element_text(size = 15, face = "bold"))
    )
} else {
  pD <- ggplot() +
    annotate("text",
      x = 0.5, y = 0.5,
      label = "fry: insufficient data for barcode panels",
      size = 5, color = "grey50"
    ) +
    theme_void()
}

save_panel(pD, file.path(RPT, "panel_D_fry_MAIN"), PD_W, PD_H)

message("F05 Panel D (fry barcode) done")
