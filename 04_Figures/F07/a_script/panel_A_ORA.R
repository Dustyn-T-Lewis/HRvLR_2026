# F07 Panel A ORA: Scatter + Flanking ORA Bars Composite
# Training concordance: HR vs LR training responses.
# X-axis = logFC_Training_HR, Y-axis = logFC_Training_LR.
# Threshold-free ORA on all proteins per quadrant, displayed as bars flanking
# the concordance scatter. Red = same-direction (concordant), Blue = opposing.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F07/a_script/style.R")
library(tidyverse)
library(fgsea)
library(ggrepel)
library(patchwork)

RPT <- "04_Figures/F07/b_reports"
DAT <- "04_Figures/F07/c_data"
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED <- unname(DIR_COLORS["Up"])
COMP_BLUE <- unname(DIR_COLORS["Down"])
N_SHOW <- 5

# -- Data --------------------------------------------------------------------
dep_df <- read_csv("03_DEP/a_non_imputed/c_data/03_combined_results.csv", show_col_types = FALSE)

# Per-protein MAR/MNAR classification no longer exists; the imputed border
# highlight falls back to the non-imputed style for all proteins.
scatter_df <- dep_df %>%
  transmute(gene,
    logFC_HR  = logFC_Training_HR,
    logFC_LR  = logFC_Training_LR,
    pi_HR     = pi_score_Training_HR,
    pi_LR     = pi_score_Training_LR,
    pi_Int    = pi_score_Training_Interaction
  ) %>%
  filter(!is.na(logFC_HR), !is.na(logFC_LR)) %>%
  mutate(
    imputed = FALSE,
    sig_class = classify_proteins_f7(pi_HR, pi_LR, pi_Int),
    is_sig = sig_class != "NS",
    quadrant = case_when(
      logFC_HR > 0 & logFC_LR > 0 ~ "Concordant Up",
      logFC_HR < 0 & logFC_LR < 0 ~ "Concordant Down",
      logFC_HR > 0 & logFC_LR < 0 ~ "Discordant (HR Up / LR Down)",
      TRUE ~ "Discordant (HR Down / LR Up)"
    )
  )

universe <- scatter_df$gene
message(sprintf(
  "  Total proteins: %d | Significant: %d",
  nrow(scatter_df), sum(scatter_df$is_sig)
))

# -- ORA ---------------------------------------------------------------------
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500)

run_set_ora <- function(genes, set_name) {
  if (length(genes) < 5) {
    return(tibble())
  }
  res <- tryCatch(
    run_ora_deduplicated(
      genes = genes, universe = universe,
      pathways = pw_collection, jaccard_cutoff = 0.5,
      min_size = 15, max_size = 500, padj_cutoff = 1
    ),
    error = function(e) {
      message("  ORA error: ", e$message)
      tibble()
    }
  )
  if (nrow(res) == 0) {
    return(tibble())
  }
  res %>%
    mutate(
      set = set_name, pathway_label = clean_pathway_name(pathway),
      neg_log10_padj = -log10(padj), significant = padj < 0.05
    ) %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = N_SHOW)
}

message("\n--- Quadrant ORA (threshold-free) ---")
ora_q1 <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Concordant Up"],
  "Concordant Up"
)
ora_q2 <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Discordant (HR Down / LR Up)"],
  "Discordant (HR Down / LR Up)"
)
ora_q3 <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Concordant Down"],
  "Concordant Down"
)
ora_q4 <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Discordant (HR Up / LR Down)"],
  "Discordant (HR Up / LR Down)"
)

all_quad_ora <- bind_rows(ora_q1, ora_q2, ora_q3, ora_q4)
if (nrow(all_quad_ora) > 0) {
  write_csv(all_quad_ora, file.path(DAT, "panel_A", "ora_quadrant.csv"))
}

# -- Scatter panel -----------------------------------------------------------
xlim_range <- c(-3.1, 3.1)
ylim_range <- c(-3.1, 3.1)

ns_df <- filter(scatter_df, sig_class == "NS")
sig_df <- filter(scatter_df, sig_class != "NS")

q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_HR > 0 & logFC_LR > 0 ~ "Q1",
    logFC_HR < 0 & logFC_LR < 0 ~ "Q3",
    logFC_HR > 0 & logFC_LR < 0 ~ "Q4",
    TRUE ~ "Q2"
  ))
q_counts <- q_df %>%
  count(q) %>%
  deframe()
q_sig <- q_df %>%
  filter(sig_class != "NS") %>%
  count(q) %>%
  deframe()
for (qq in c("Q1", "Q2", "Q3", "Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- sig_df %>%
  group_by(sig_class) %>%
  arrange(desc(abs(logFC_HR) + abs(logFC_LR))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    label_fill = SIG_LABEL_FILL_F7[as.character(sig_class)],
    label_text_col = SIG_LABEL_TEXT_F7[as.character(sig_class)]
  )

txt_gene <- scale_text(BASE_GENE, 200) * 1.2
txt_quad <- scale_text(BASE_QUADRANT, 200) * 1.3

# Center-axis tick labels
x_breaks <- seq(-3, 3, 1)
y_breaks <- seq(-3, 3, 1)
x_tick_df <- tibble(
  x = x_breaks[x_breaks != 0], y = 0,
  label = as.character(x_breaks[x_breaks != 0])
)
y_tick_df <- tibble(
  x = 0, y = y_breaks[y_breaks != 0],
  label = as.character(y_breaks[y_breaks != 0])
)

p_scatter <- ggplot(mapping = aes(x = logFC_HR, y = logFC_LR)) +
  annotate("rect",
    xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
    fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
    fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
    fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
    fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "black", linewidth = 0.3
  ) +
  geom_text(
    data = x_tick_df, aes(x = x, y = y, label = label),
    vjust = 1.5, size = 3.2, color = "grey40", fontface = "bold"
  ) +
  geom_text(
    data = y_tick_df, aes(x = x, y = y, label = label),
    hjust = 1.5, size = 3.2, color = "grey40", fontface = "bold"
  ) +
  geom_point(
    data = ns_df, color = "grey80", fill = "grey85", shape = 21,
    size = 1.0, alpha = 0.3, stroke = 0.2
  ) +
  geom_point(
    data = sig_df, aes(fill = sig_class), shape = 21,
    size = ifelse(sig_df$sig_class == "NS", 1.8, 2.3),
    color = ifelse(sig_df$imputed, "black", "grey75"),
    alpha = case_when(
      sig_df$sig_class == "NS" ~ 0.30,
      sig_df$sig_class == "Interaction" ~ 0.55,
      sig_df$sig_class == "Sig Both" ~ 0.75,
      TRUE ~ 0.85
    ),
    stroke = ifelse(sig_df$sig_class == "NS", 0.6, 0.9)
  ) +
  scale_fill_manual(values = SIG_COLORS_F7, name = "Significance") +
  geom_label_repel(
    data = label_df, aes(label = gene),
    fill = label_df$label_fill, color = label_df$label_text_col,
    size = txt_gene, fontface = "italic", max.overlaps = 50,
    segment.size = 0.2, segment.color = "grey50",
    min.segment.length = 0, show.legend = FALSE,
    box.padding = 0.4, point.padding = 0.3,
    force = 4, force_pull = 0.3,
    label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
    label.size = 0.15, seed = 42,
    xlim = c(-3, 3) * 0.92, ylim = c(-3, 3) * 0.92
  ) +
  # Quadrant labels
  annotate("label",
    x = xlim_range[2], y = ylim_range[2],
    label = sprintf("Concordant Up  %s/%s", q_sig["Q1"], q_counts["Q1"]),
    hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
    color = COMP_RED, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt")
  ) +
  annotate("label",
    x = xlim_range[1], y = ylim_range[1],
    label = sprintf("Concordant Down  %s/%s", q_sig["Q3"], q_counts["Q3"]),
    hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
    color = COMP_RED, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt")
  ) +
  annotate("label",
    x = xlim_range[1], y = ylim_range[2],
    label = sprintf("Discordant (HR\u2193 LR\u2191)  %s/%s", q_sig["Q2"], q_counts["Q2"]),
    hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
    color = COMP_BLUE, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt")
  ) +
  annotate("label",
    x = xlim_range[2], y = ylim_range[1],
    label = sprintf("Discordant (HR\u2191 LR\u2193)  %s/%s", q_sig["Q4"], q_counts["Q4"]),
    hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
    color = COMP_BLUE, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt")
  ) +
  # Axis titles
  annotate("text",
    x = 2.5, y = 0,
    label = expression(log[2] * FC ~ "(Training HR)"),
    hjust = 0.5, vjust = -0.4, size = 3.2, color = "grey30", fontface = "bold"
  ) +
  annotate("text",
    x = 0, y = 2.5,
    label = expression(log[2] * FC ~ "(Training LR)"),
    hjust = 0.5, vjust = -0.4, size = 3.2, color = "grey30", fontface = "bold",
    angle = 90
  ) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(2, 0, 2, 0, "mm"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(5, "mm"),
    legend.margin = margin(-5, 0, 0, 0),
    legend.box.margin = margin(-8, 0, 0, 0)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 5, alpha = 0.8)))

# -- Half-bar builder --------------------------------------------------------
make_half_bars <- function(df, fill_color, side, ylim,
                           display_labels = character(0)) {
  bar_h <- 0.24
  n_bars <- if (is.null(df) || nrow(df) == 0) 0L else min(nrow(df), 5L)

  if (n_bars == 0) {
    return(ggplot() +
      theme_void() +
      scale_y_continuous(limits = ylim, expand = c(0, 0)))
  }

  y_pos <- if (ylim[1] >= 0) {
    rev(seq(0.3, 2.7, length.out = 5))[seq_len(n_bars)]
  } else {
    seq(-0.3, -2.7, length.out = 5)[seq_len(n_bars)]
  }

  bars <- df %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = 5) %>%
    mutate(
      y = y_pos,
      bar_fill = ifelse(significant, scales::alpha(fill_color, 0.85),
        scales::alpha(fill_color, 0.30)
      ),
      display_name = ifelse(pathway_label %in% names(display_labels),
        display_labels[pathway_label], pathway_label
      ),
      display_name = ifelse(!grepl("\n", display_name),
        stringr::str_wrap(display_name, width = 35),
        display_name
      ),
      star_raw = sig_stars(padj),
      star = gsub("\\*", "*\n", star_raw) %>% sub("\n$", "", .)
    )

  x_max <- max(bars$neg_log10_padj)
  x_display_max <- x_max * 1.08

  bars <- bars %>%
    mutate(
      display_name = ifelse(!significant & nchar(display_name) > 18,
        stringr::str_wrap(display_name, width = 14),
        display_name
      ),
      text_size = case_when(
        significant & nchar(display_name) > 28 & (neg_log10_padj / x_max) < 0.3 ~ 3.0,
        significant & nchar(display_name) > 28 ~ 3.8,
        significant ~ 4.5,
        !significant & nchar(display_name) > 18 ~ 2.6,
        !significant ~ 3.0,
        TRUE ~ 3.2
      )
    )

  is_upper <- ylim[1] >= 0
  brk_fn <- function(limits) {
    b <- scales::pretty_breaks(n = 3)(limits)
    b[b != 0]
  }

  p <- ggplot(bars, aes(y = y)) +
    geom_rect(
      aes(
        xmin = 0, xmax = neg_log10_padj,
        ymin = y - bar_h / 2, ymax = y + bar_h / 2
      ),
      fill = bars$bar_fill, color = "black", linewidth = 0.3
    ) +
    geom_text(
      data = bars %>% filter(significant),
      aes(x = neg_log10_padj / 2, y = y, label = display_name),
      hjust = 0.5, size = (bars %>% filter(significant))$text_size,
      fontface = "bold", color = "white", lineheight = 0.85
    ) +
    geom_text(
      data = bars %>% filter(!significant),
      aes(x = neg_log10_padj / 2, y = y, label = display_name),
      hjust = 0.5, size = (bars %>% filter(!significant))$text_size,
      fontface = "bold", color = "grey30", lineheight = 0.85
    ) +
    geom_text(aes(x = neg_log10_padj + x_max * 0.025, label = star),
      hjust = 0.5, vjust = 0.5,
      size = 3.5, fontface = "bold", color = "black",
      lineheight = 0.7
    ) +
    annotate("segment",
      x = 0, xend = x_display_max, y = -Inf, yend = -Inf,
      color = "grey40", linewidth = 0.3
    ) +
    labs(
      x = if (!is_upper) expression(-log[10](p[adj])) else NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 8.5, face = "bold"),
      axis.title.x = if (!is_upper) {
        element_text(size = 9.5, face = "bold")
      } else {
        element_blank()
      },
      axis.line.x = element_blank(),
      axis.ticks.x = element_line(color = "grey40", linewidth = 0.3),
      plot.margin = if (side == "left") {
        margin(2, 0, 2, 3, "mm")
      } else {
        margin(2, 3, 2, 0, "mm")
      }
    )

  if (side == "left") {
    p + scale_x_reverse(
      limits = c(x_display_max, 0),
      breaks = brk_fn,
      expand = expansion(mult = c(0.02, 0))
    ) +
      scale_y_continuous(limits = ylim, expand = c(0, 0))
  } else {
    p + scale_x_continuous(
      limits = c(0, x_display_max),
      breaks = brk_fn,
      expand = expansion(mult = c(0, 0.02))
    ) +
      scale_y_continuous(limits = ylim, expand = c(0, 0))
  }
}

# 4 half-bar panels
p_ul <- make_half_bars(ora_q2, scales::alpha(COMP_BLUE, 0.30), "left", c(0, 3.1),
  display_labels = DISPLAY_LABELS_F07
)
p_ll <- make_half_bars(ora_q3, scales::alpha(COMP_RED, 0.30), "left", c(-3.1, 0),
  display_labels = DISPLAY_LABELS_F07
)
p_ur <- make_half_bars(ora_q1, scales::alpha(COMP_RED, 0.30), "right", c(0, 3.1),
  display_labels = DISPLAY_LABELS_F07
)
p_lr <- make_half_bars(ora_q4, scales::alpha(COMP_BLUE, 0.30), "right", c(-3.1, 0),
  display_labels = DISPLAY_LABELS_F07
)

# -- Composite ---------------------------------------------------------------
design <- c(area(1, 1), area(1, 2, 2, 2), area(1, 3), area(2, 1), area(2, 3))
n_total <- nrow(scatter_df)
n_sig <- sum(scatter_df$is_sig)
n_enrich <- if (nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_pear <- cor(scatter_df$logFC_HR, scatter_df$logFC_LR, use = "complete.obs")

composite <- p_ul + p_scatter + p_ur + p_ll + p_lr +
  plot_layout(design = design, widths = c(1.4, 2, 1.4)) +
  plot_annotation(
    title = "Training Concordance: Quadrant ORA",
    subtitle = sprintf(
      "Threshold-free ORA (hypergeometric) | N = %d | %d DEPs (Pi < 0.05) | %d enriched (FDR < 0.05) | r = %.2f",
      n_total, n_sig, n_enrich, r_pear
    ),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey30")
    )
  )

COMP_W <- 450
COMP_H <- 220
ggsave(file.path(RPT, "panel_A_ORA_composite_MAIN.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_device
)
ggsave(file.path(RPT, "panel_A_ORA_composite_MAIN.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

message("\nF07 Panel A ORA composite done")
