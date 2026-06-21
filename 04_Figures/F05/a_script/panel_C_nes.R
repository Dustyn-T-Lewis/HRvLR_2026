# F05 Panel C: NES Scatter (GO Slim + Hallmark, a priori)
# Two side-by-side: Training_HR NES vs Training_LR NES, Acute_HR NES vs Acute_LR NES
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F05/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

PC_W <- 400; PC_H <- 200

RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# -- Load fGSEA cache ---------------------------------------------------------

fgsea_cache <- file.path(DAT, "fgsea_cache.csv")
if (!file.exists(fgsea_cache)) stop("fGSEA cache not found -- run prepare_data.R first")
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

# -- Pair definitions ----------------------------------------------------------

pairs <- list(
  Training = list(
    ctr_x = "Training_HR", ctr_y = "Training_LR",
    label_x = "HR", label_y = "LR", title = "Training"
  ),
  Acute = list(
    ctr_x = "Acute_HR", ctr_y = "Acute_LR",
    label_x = "HR", label_y = "LR", title = "Acute"
  )
)

# -- Significance colors (shared across panels) --------------------------------

SIG_COLORS <- c(
  "Sig Both"    = "#2E7D32",
  "Sig HR only" = "#4393C3",
  "Sig LR only" = "#D6604D",
  "NS"          = "grey70"
)

SIG_LABEL_FILL <- c(
  "Sig Both"    = alpha("#2E7D32", 0.15),
  "Sig HR only" = alpha("#4393C3", 0.15),
  "Sig LR only" = alpha("#D6604D", 0.15),
  "NS"          = alpha("grey70", 0.15)
)

SIG_LABEL_TEXT <- c(
  "Sig Both"    = "#2E7D32",
  "Sig HR only" = "#4393C3",
  "Sig LR only" = "#D6604D",
  "NS"          = "grey50"
)

# -- Build one NES scatter panel -----------------------------------------------

make_nes_scatter <- function(fgsea, pair, half_w) {
  # Filter to GO Slim + Hallmark (a priori)
  # Use "GO:BP" as proxy for GO Slim when full GO Slim label not available
  fgsea_hg <- fgsea %>%
    filter(database %in% c("Hallmark", "GO:BP"),
           contrast %in% c(pair$ctr_x, pair$ctr_y))

  # Pivot to wide
  fgsea_wide <- fgsea_hg %>%
    select(pathway, contrast, NES, padj, size, database) %>%
    pivot_wider(id_cols = c(pathway, database), names_from = contrast,
                values_from = c(NES, padj, size)) %>%
    filter(!is.na(.data[[paste0("NES_", pair$ctr_x)]]),
           !is.na(.data[[paste0("NES_", pair$ctr_y)]]))

  nes_x_col  <- paste0("NES_", pair$ctr_x)
  nes_y_col  <- paste0("NES_", pair$ctr_y)
  padj_x_col <- paste0("padj_", pair$ctr_x)
  padj_y_col <- paste0("padj_", pair$ctr_y)
  size_col   <- paste0("size_", pair$ctr_x)

  fgsea_wide <- fgsea_wide %>%
    mutate(
      nes_x = .data[[nes_x_col]],
      nes_y = .data[[nes_y_col]],
      padj_x = .data[[padj_x_col]],
      padj_y = .data[[padj_y_col]],
      set_size = coalesce(.data[[size_col]],
                           .data[[paste0("size_", pair$ctr_y)]]),
      sig_x = !is.na(padj_x) & padj_x < 0.05,
      sig_y = !is.na(padj_y) & padj_y < 0.05,
      significance = case_when(
        sig_x & sig_y ~ "Sig Both",
        sig_x         ~ paste("Sig", pair$label_x, "only"),
        sig_y         ~ paste("Sig", pair$label_y, "only"),
        TRUE          ~ "NS"
      ) %>% factor(levels = names(SIG_COLORS)),
      pathway_label = clean_pathway_name(pathway)
    )

  fgsea_sig <- fgsea_wide %>% filter(significance != "NS")
  n_total   <- nrow(fgsea_wide)
  n_sig     <- nrow(fgsea_sig)

  # Spearman on all
  nes_cor_all <- cor.test(fgsea_wide$nes_x, fgsea_wide$nes_y, method = "spearman")
  nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, n_total)

  # Spearman on sig
  nes_cor_sig <- if (n_sig >= 3) {
    cor.test(fgsea_sig$nes_x, fgsea_sig$nes_y, method = "spearman")
  } else NULL

  # Pathway concordance
  n_conc <- sum(fgsea_sig$nes_x > 0 & fgsea_sig$nes_y > 0) +
            sum(fgsea_sig$nes_x < 0 & fgsea_sig$nes_y < 0)
  pw_conc_frac <- if (n_sig > 0) n_conc / n_sig else 0
  pw_conc_binom <- if (n_sig > 0) binom.test(n_conc, n_sig) else NULL
  pw_conc_ci <- if (!is.null(pw_conc_binom)) pw_conc_binom$conf.int * 100 else c(NA, NA)

  # Quadrant counts (sig only)
  n_conc_tr <- sum(fgsea_sig$nes_x > 0 & fgsea_sig$nes_y > 0)
  n_conc_bl <- sum(fgsea_sig$nes_x < 0 & fgsea_sig$nes_y < 0)
  n_disc_q2 <- sum(fgsea_sig$nes_x < 0 & fgsea_sig$nes_y > 0)
  n_disc_q4 <- sum(fgsea_sig$nes_x > 0 & fgsea_sig$nes_y < 0)

  nes_lim <- max(abs(c(fgsea_wide$nes_x, fgsea_wide$nes_y))) * 1.15

  # Sig labels on all sig terms
  label_pw <- fgsea_sig %>%
    mutate(
      label_fill     = SIG_LABEL_FILL[as.character(significance)],
      label_text_col = SIG_LABEL_TEXT[as.character(significance)]
    ) %>%
    arrange(significance)

  ns_df <- fgsea_wide %>% filter(significance == "NS")
  sig_df <- fgsea_wide %>% filter(significance != "NS") %>%
    mutate(
      border_col = ifelse(database == "Hallmark", "black", "grey75"),
      bubble_alpha = case_when(
        significance == "Sig Both" ~ 0.75,
        TRUE ~ 0.85
      )
    )

  txt_gene <- scale_text(BASE_GENE, half_w)
  txt_quad <- scale_text(BASE_QUADRANT, half_w)

  rho_sig_str <- if (!is.null(nes_cor_sig)) {
    sprintf(", rho(sig)=%.2f", nes_cor_sig$estimate)
  } else ""

  subtitle_str <- sprintf(
    "GO:BP + Hallmark | %d pathways (%d sig) | rho(all)=%.2f [%.2f,%.2f]%s | %.0f%% conc.",
    n_total, n_sig,
    nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
    rho_sig_str, pw_conc_frac * 100)

  p <- ggplot(mapping = aes(x = nes_x, y = nes_y)) +
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
    geom_point(data = ns_df, aes(x = nes_x, y = nes_y),
               size = 1.5, fill = "grey70",
               shape = 21, color = "grey55", alpha = 0.40, stroke = 0.4) +
    geom_point(data = sig_df, aes(x = nes_x, y = nes_y,
               fill = significance, size = set_size),
               shape = 21, color = sig_df$border_col,
               alpha = sig_df$bubble_alpha, stroke = 0.8) +
    scale_fill_manual(values = SIG_COLORS, name = "Significance") +
    scale_size_continuous(range = c(2, 8), name = "Set size",
                          breaks = c(20, 50, 100, 200)) +
    geom_label_repel(data = label_pw,
                     aes(x = nes_x, y = nes_y, label = pathway_label),
                     fill = label_pw$label_fill, color = label_pw$label_text_col,
                     size = txt_gene, fontface = "bold",
                     max.overlaps = 40,
                     segment.size = 0.2, segment.color = "grey50",
                     min.segment.length = 0, show.legend = FALSE,
                     box.padding = 0.5, force = 3, force_pull = 0.5,
                     label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
                     label.size = 0.15, seed = 42) +
    annotate("label", x = Inf, y = Inf,
             label = sprintf("Concordant Up  n=%d", n_conc_tr),
             hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
             color = "#D6604D", fill = alpha("white", 0.9),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = -Inf, y = -Inf,
             label = sprintf("Concordant Down  n=%d", n_conc_bl),
             hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
             color = "#D6604D", fill = alpha("white", 0.9),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = -Inf, y = Inf,
             label = sprintf("Discordant  n=%d", n_disc_q2),
             hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
             color = "#4393C3", fill = alpha("white", 0.9),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = Inf, y = -Inf,
             label = sprintf("Discordant  n=%d", n_disc_q4),
             hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
             color = "#4393C3", fill = alpha("white", 0.9),
             label.padding = unit(2.5, "pt")) +
    scale_x_continuous(expand = expansion(0, 0)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
    labs(
      title = pair$title,
      subtitle = subtitle_str,
      x = sprintf("NES (%s)", pair$label_x),
      y = sprintf("NES (%s)", pair$label_y)
    ) +
    FIG_THEME +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 8, face = "bold"),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(3, "mm"),
      legend.margin   = margin(0, 0, 0, 0)
    ) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 0.8)),
           size = guide_legend(nrow = 1))

  # Export data
  export_df <- fgsea_wide %>%
    transmute(
      pathway, pathway_label, database,
      nes_x = round(nes_x, 3), nes_y = round(nes_y, 3),
      padj_x = signif(padj_x, 4), padj_y = signif(padj_y, 4),
      significance = as.character(significance), set_size
    ) %>%
    arrange(significance, desc(abs(nes_x) + abs(nes_y)))

  list(plot = p, export = export_df)
}

# -- Build two panels ----------------------------------------------------------

half_w <- PC_W / 2
res_train <- make_nes_scatter(fgsea_all, pairs$Training, half_w)
res_acute <- make_nes_scatter(fgsea_all, pairs$Acute,    half_w)

pC <- (res_train$plot | res_acute$plot) +
  plot_annotation(
    title = "Pathway-Level Concordance (fGSEA NES)",
    tag_levels = list(c("C", "")),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.tag   = element_text(size = 15, face = "bold")
    )
  )

save_panel(pC, file.path(RPT, "panel_C_nes_MAIN"), PC_W, PC_H)

# Export CSVs
write_csv(res_train$export, file.path(DAT, "nes_scatter_training.csv"))
write_csv(res_acute$export,  file.path(DAT, "nes_scatter_acute.csv"))

message("F05 Panel C (NES scatter) done")
