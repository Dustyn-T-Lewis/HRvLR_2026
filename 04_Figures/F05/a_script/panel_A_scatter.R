# F05 Panel A: Between-Group Concordance Scatters (Training + Acute)
# Two side-by-side scatters: Training_HR vs Training_LR, Acute_HR vs Acute_LR
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F05/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(boot)
})

PA_W <- 400; PA_H <- 200

RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

# -- Pair definitions ----------------------------------------------------------

pairs <- list(
  Training = list(
    lfc_x = "logFC_Training_HR", lfc_y = "logFC_Training_LR",
    pi_x  = "pi_score_Training_HR", pi_y = "pi_score_Training_LR",
    t_x   = "t_Training_HR",  t_y   = "t_Training_LR",
    label_x = "HR", label_y = "LR",
    title = "Training"
  ),
  Acute = list(
    lfc_x = "logFC_Acute_HR", lfc_y = "logFC_Acute_LR",
    pi_x  = "pi_score_Acute_HR", pi_y = "pi_score_Acute_LR",
    t_x   = "t_Acute_HR", t_y   = "t_Acute_LR",
    label_x = "HR", label_y = "LR",
    title = "Acute"
  )
)

# -- Build one scatter panel ---------------------------------------------------

make_scatter <- function(dep, pair, half_w) {
  sdf <- dep %>%
    transmute(
      gene,
      lfc_x = .data[[pair$lfc_x]], lfc_y = .data[[pair$lfc_y]],
      pi_x  = .data[[pair$pi_x]],  pi_y  = .data[[pair$pi_y]]
    ) %>%
    filter(!is.na(lfc_x), !is.na(lfc_y)) %>%
    mutate(
      significance = classify_concordance(pi_x, pi_y, pair$label_x, pair$label_y),
      quadrant = classify_quadrant(lfc_x, lfc_y, pair$label_x, pair$label_y)
    )

  conc_cols <- make_concordance_colors(pair$label_x, pair$label_y)
  n_obs <- nrow(sdf)

  # Correlations
  cor_r   <- cor.test(sdf$lfc_x, sdf$lfc_y, method = "pearson",  conf.level = 0.95)
  cor_rho <- cor.test(sdf$lfc_x, sdf$lfc_y, method = "spearman", conf.level = 0.95)
  rho_ci  <- fisher_z_ci(cor_rho$estimate, n_obs)
  sign_conc <- mean(sign(sdf$lfc_x) == sign(sdf$lfc_y)) * 100

  # Sign concordance bootstrap CI
  set.seed(42)
  boot_sc <- boot::boot(sdf, function(d, i)
    mean(sign(d$lfc_x[i]) == sign(d$lfc_y[i])) * 100, R = 5000)
  sc_ci <- tryCatch(
    boot::boot.ci(boot_sc, type = "bca", conf = 0.95)$bca[4:5],
    error = function(e) quantile(boot_sc$t, c(0.025, 0.975)))

  # Sig-only stats
  sig_mask <- sdf$significance != "NS"
  n_sig    <- sum(sig_mask)
  conc_sig <- if (n_sig > 3) mean(sign(sdf$lfc_x[sig_mask]) == sign(sdf$lfc_y[sig_mask])) * 100 else NA_real_
  rho_sig  <- if (n_sig > 3) cor(sdf$lfc_x[sig_mask], sdf$lfc_y[sig_mask], method = "spearman") else NA_real_
  rho_sig_ci <- if (n_sig > 3) fisher_z_ci(rho_sig, n_sig) else c(lo = NA, hi = NA)

  # Quadrant counts
  q_df <- sdf %>%
    mutate(q = case_when(
      lfc_x > 0 & lfc_y > 0 ~ "Q1",
      lfc_x < 0 & lfc_y < 0 ~ "Q3",
      lfc_x > 0 & lfc_y < 0 ~ "Q4",
      TRUE ~ "Q2"
    ))
  q_counts <- q_df %>% count(q) %>% deframe()
  q_sig    <- q_df %>% filter(significance != "NS") %>% count(q) %>% deframe()
  for (qq in c("Q1", "Q2", "Q3", "Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

  # Axis limits
  max_abs <- max(abs(c(sdf$lfc_x, sdf$lfc_y)), na.rm = TRUE) * 1.1
  lim <- c(-max_abs, max_abs)

  # Label top proteins
  label_df <- sdf %>%
    filter(significance != "NS") %>%
    group_by(significance) %>%
    arrange(desc(abs(lfc_x) + abs(lfc_y))) %>%
    slice_head(n = 4) %>%
    ungroup()

  ns_df  <- sdf %>% filter(significance == "NS")
  sig_df <- sdf %>% filter(significance != "NS")

  txt_gene <- scale_text(BASE_GENE, half_w)
  txt_quad <- scale_text(BASE_QUADRANT, half_w)
  txt_stat <- scale_text(BASE_STAT, half_w)

  # Subtitle
  sub_str <- sprintf(
    "All (n=%s): \u03c1=%.2f [%.2f,%.2f], r=%.2f [%.2f,%.2f] | conc.=%.0f%%\nSig (n=%d): \u03c1=%.2f [%.2f,%.2f] | conc.=%.0f%%",
    format(n_obs, big.mark = ","),
    cor_rho$estimate, rho_ci[1], rho_ci[2],
    cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2],
    sign_conc,
    n_sig,
    ifelse(is.na(rho_sig), 0, rho_sig),
    ifelse(is.na(rho_sig_ci["lo"]), 0, rho_sig_ci["lo"]),
    ifelse(is.na(rho_sig_ci["hi"]), 0, rho_sig_ci["hi"]),
    ifelse(is.na(conc_sig), 0, conc_sig))

  p <- ggplot(mapping = aes(x = lfc_x, y = lfc_y)) +
    # Concordant quadrants (pink)
    annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
             fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
    # Discordant quadrants (blue)
    annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
             fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
             fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.3) +
    geom_point(data = ns_df, size = 1.0, fill = "grey85", shape = 21,
               color = "grey80", alpha = 0.30, stroke = 0.2) +
    geom_point(data = sig_df, aes(fill = significance), shape = 21,
               size = 2.0, color = "grey75", alpha = 0.75, stroke = 0.7) +
    scale_fill_manual(values = conc_cols, name = "Significance") +
    geom_label_repel(data = label_df, aes(label = gene),
                     size = txt_gene, fontface = "italic",
                     fill = alpha("white", 0.85), color = "grey20",
                     max.overlaps = 30,
                     segment.size = 0.2, segment.color = "grey50",
                     min.segment.length = 0, show.legend = FALSE,
                     box.padding = 0.5, force = 3, force_pull = 0.5,
                     label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
                     label.size = 0.15, seed = 42) +
    # Quadrant counts
    annotate("label", x = Inf, y = Inf,
             label = sprintf("Concordant Up\u2002n=%s/%s", q_sig["Q1"], q_counts["Q1"]),
             hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
             color = "#D6604D", fill = alpha("white", 0.9),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = -Inf, y = -Inf,
             label = sprintf("Concordant Down\u2002n=%s/%s", q_sig["Q3"], q_counts["Q3"]),
             hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
             color = "#D6604D", fill = alpha("white", 0.9),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = -Inf, y = Inf,
             label = sprintf("Discordant\u2002n=%s/%s", q_sig["Q2"], q_counts["Q2"]),
             hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
             color = "#4393C3", fill = alpha("white", 0.9),
             label.padding = unit(2.5, "pt")) +
    annotate("label", x = Inf, y = -Inf,
             label = sprintf("Discordant\u2002n=%s/%s", q_sig["Q4"], q_counts["Q4"]),
             hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
             color = "#4393C3", fill = alpha("white", 0.9),
             label.padding = unit(2.5, "pt")) +
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, expand = FALSE) +
    labs(
      title = pair$title,
      subtitle = sub_str,
      x = bquote(log[2]*FC ~ "(" * .(pair$label_x) * ")"),
      y = bquote(log[2]*FC ~ "(" * .(pair$label_y) * ")")
    ) +
    FIG_THEME +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 8, face = "bold"),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(3, "mm"),
      legend.margin   = margin(0, 0, 0, 0)
    ) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 0.8)))

  p
}

# -- Build two panels ----------------------------------------------------------

half_w <- PA_W / 2
p_train <- make_scatter(dep_df, pairs$Training, half_w)
p_acute <- make_scatter(dep_df, pairs$Acute,    half_w)

pA <- (p_train | p_acute) +
  plot_annotation(
    title = "Protein-Level Concordance: HR vs LR",
    subtitle = "Non-imputed limma + duplicateCorrelation | Pi-score significance (< 0.05)",
    tag_levels = list(c("A", "")),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, face = "bold.italic", color = "grey30"),
      plot.tag      = element_text(size = 15, face = "bold")
    )
  )

save_panel(pA, file.path(RPT, "panel_A_scatter_MAIN"), PA_W, PA_H)

message("F05 Panel A (scatter) done")
