# F06 Panel B: RRHO2 Threshold-Free Concordance Heatmaps
# Two side-by-side RRHO2 maps: within HR (training vs acute) + within LR (training vs acute)
# Stratified hypergeometric via RRHO2 package (Cahill et al. 2018), stepsize = 20
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(RRHO2)
})

PB_W <- 400; PB_H <- 220

RPT <- "04_Figures/F06/b_reports"
DAT <- "04_Figures/F06/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

JET_COLORS <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                "yellow", "#FF7F00", "red", "#7F0000")

# -- Pair definitions ----------------------------------------------------------

pairs <- list(
  Training = list(
    t_x = "t_Training_HR", t_y = "t_Acute_HR",
    label_x = "Training", label_y = "Acute", title = "High Responders"
  ),
  Acute = list(
    t_x = "t_Training_LR", t_y = "t_Acute_LR",
    label_x = "Training", label_y = "Acute", title = "Low Responders"
  )
)

# -- Build one RRHO2 heatmap panel ---------------------------------------------

make_rrho2 <- function(dep, pair, half_w) {
  rr_df <- dep %>%
    transmute(gene, t_x = .data[[pair$t_x]], t_y = .data[[pair$t_y]]) %>%
    filter(!is.na(t_x) & !is.na(t_y)) %>%
    distinct(gene, .keep_all = TRUE)

  n_shared <- nrow(rr_df)

  list1 <- data.frame(gene = rr_df$gene, score = rr_df$t_x, stringsAsFactors = FALSE)
  list2 <- data.frame(gene = rr_df$gene, score = rr_df$t_y, stringsAsFactors = FALSE)

  rrho_obj <- RRHO2_initialize(
    list1, list2,
    labels          = c(pair$label_x, pair$label_y),
    log10.ind       = TRUE,
    multipleTesting = "none",
    boundary        = 0.02,
    method          = "hyper",
    stepsize        = 20
  )

  hmat <- rrho_obj$hypermat
  nr <- nrow(hmat); nc <- ncol(hmat)

  na_rows <- which(apply(hmat, 1, function(r) all(is.na(r))))
  na_cols <- which(apply(hmat, 2, function(c) all(is.na(c))))

  if (length(na_rows) > 0 && length(na_cols) > 0) {
    row_before <- 1:(min(na_rows) - 1)
    row_after  <- (max(na_rows) + 1):nr
    col_before <- 1:(min(na_cols) - 1)
    col_after  <- (max(na_cols) + 1):nc
  } else {
    mid <- floor(nr / 2)
    row_before <- 1:mid
    row_after  <- (mid + 1):nr
    col_before <- 1:mid
    col_after  <- (mid + 1):nc
  }

  max_UU <- max(hmat[row_before, col_before], na.rm = TRUE)
  max_DD <- max(hmat[row_after,  col_after],  na.rm = TRUE)
  max_UD <- max(hmat[row_before, col_after],  na.rm = TRUE)
  max_DU <- max(hmat[row_after,  col_before], na.rm = TRUE)

  hotspot_UU <- rrho_obj$genelist_uu$gene_list_overlap_uu
  hotspot_DD <- rrho_obj$genelist_dd$gene_list_overlap_dd
  hotspot_UD <- rrho_obj$genelist_ud$gene_list_overlap_ud
  hotspot_DU <- rrho_obj$genelist_du$gene_list_overlap_du

  n_UU <- length(hotspot_UU)
  n_DD <- length(hotspot_DD)
  n_UD <- length(hotspot_UD)
  n_DU <- length(hotspot_DU)

  hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
    mutate(neg_log10_p = as.vector(hmat))

  max_val <- max(hmat_df$neg_log10_p, na.rm = TRUE)

  txt_quad <- scale_text(BASE_QUADRANT, half_w)
  txt_stat <- scale_text(BASE_STAT, half_w) * 0.75

  ann_x_left  <- mean(row_before)
  ann_x_right <- mean(row_after)
  ann_y_bot   <- mean(col_before)
  ann_y_top   <- mean(col_after)

  p <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_p)) +
    geom_raster() +
    scale_fill_gradientn(
      colors   = JET_COLORS,
      limits   = c(0, max_val),
      na.value = "white",
      name     = expression(-log[10](P)),
      guide    = guide_colorbar(
        barwidth = unit(20, "mm"), barheight = unit(2.5, "mm"),
        title.position = "left", title.hjust = 1,
        title.theme = element_text(size = 7, face = "bold"))) +
    annotate("text", x = ann_x_left, y = ann_y_bot,
             label = "Concordant Up",
             color = "white", fontface = "bold", size = txt_quad) +
    annotate("text", x = ann_x_right, y = ann_y_top,
             label = "Concordant Down",
             color = "white", fontface = "bold", size = txt_quad) +
    annotate("text", x = ann_x_left, y = ann_y_top,
             label = paste0("Discordant\n", pair$label_x, " Up / ", pair$label_y, " Down"),
             color = "white", fontface = "bold", size = txt_quad) +
    annotate("text", x = ann_x_right, y = ann_y_bot,
             label = paste0("Discordant\n", pair$label_x, " Down / ", pair$label_y, " Up"),
             color = "white", fontface = "bold", size = txt_quad) +
    annotate("text", x = ann_x_left, y = ann_y_bot - diff(range(col_before)) * 0.18,
             label = sprintf("max=%.1f | n=%d", max_UU, n_UU),
             color = "white", size = txt_stat) +
    annotate("text", x = ann_x_right, y = ann_y_top - diff(range(col_after)) * 0.15,
             label = sprintf("max=%.1f | n=%d", max_DD, n_DD),
             color = "white", size = txt_stat) +
    annotate("text", x = ann_x_left, y = ann_y_top - diff(range(col_after)) * 0.18,
             label = sprintf("max=%.1f | n=%d", max_UD, n_UD),
             color = "white", size = txt_stat) +
    annotate("text", x = ann_x_right, y = ann_y_bot - diff(range(col_before)) * 0.18,
             label = sprintf("max=%.1f | n=%d", max_DU, n_DU),
             color = "white", size = txt_stat) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = pair$title,
      subtitle = sprintf("Stratified hypergeometric | %d genes | stepsize=20",
                          n_shared),
      x = bquote(.(pair$label_x) ~ "rank" ~ (Up %->% Down)),
      y = bquote(.(pair$label_y) ~ "rank" ~ (Up %->% Down))
    ) +
    FIG_THEME +
    theme(
      axis.text        = element_blank(),
      axis.title.x     = element_text(margin = margin(t = 2)),
      axis.title.y     = element_text(margin = margin(r = 2)),
      axis.ticks       = element_blank(),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      legend.position  = "bottom",
      legend.margin    = margin(0, 0, 0, 0),
      plot.margin = margin(2, 2, 2, 2, "mm")
    ) +
    coord_fixed(ratio = 1)

  meta <- tibble(
    pair = pair$title,
    max_concordant_up   = round(max_UU, 2),
    max_concordant_down = round(max_DD, 2),
    max_discordant_1    = round(max_UD, 2),
    max_discordant_2    = round(max_DU, 2),
    n_hotspot_UU = n_UU, n_hotspot_DD = n_DD,
    n_hotspot_UD = n_UD, n_hotspot_DU = n_DU,
    n_shared = n_shared, stepsize = 20, grid_size = nr
  )

  list(plot = p, meta = meta)
}

# -- Build two panels ----------------------------------------------------------

half_w <- PB_W / 2
res_train <- make_rrho2(dep_df, pairs$Training, half_w)
res_acute <- make_rrho2(dep_df, pairs$Acute,    half_w)

pB <- (res_train$plot | res_acute$plot) +
  plot_annotation(
    title = "Threshold-Free Concordance (RRHO2)",
    tag_levels = list(c("B", "")),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.tag   = element_text(size = 15, face = "bold")
    )
  )

save_panel(pB, file.path(RPT, "panel_B_RRHO2_MAIN"), PB_W, PB_H)

rrho2_meta <- bind_rows(res_train$meta, res_acute$meta)
write_csv(rrho2_meta, file.path(DAT, "rrho2_summary.csv"))

message("F06 Panel B (RRHO2) done")
