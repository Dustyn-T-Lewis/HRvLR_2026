# Figure 1 — Panel B: Effect Size Distribution (logFC Histograms) ----
# Layout: Baseline (top, full-width) | Training HR+LR | Acute HR+LR

if (!exists("meta")) source("04_Figures/F1/a_script/HRvLR_F1_setup.R")

PB_W <- 130
PB_H <- 120

# ── Reshape to long format (signed logFC, clipped to [-1, 1]) ----

lfc_long <- dep_df |>
  dplyr::select(gene, starts_with("logFC_")) |>
  tidyr::pivot_longer(starts_with("logFC_"),
                      names_to = "contrast", values_to = "logFC") |>
  dplyr::mutate(contrast = stringr::str_remove(contrast, "logFC_")) |>
  dplyr::filter(!is.na(logFC), contrast %in% MAIN_CONTRASTS) |>
  dplyr::mutate(abs_logFC = abs(logFC),
                contrast  = factor(contrast, levels = MAIN_CONTRASTS))

# ── Summary stats per contrast ----
set.seed(42)

lfc_stats <- lfc_long |>
  dplyr::group_by(contrast) |>
  dplyr::summarise(
    med_lfc     = median(logFC),
    med_abs_lfc = median(abs_logFC),
    ci          = list(boot_median_ci(abs_logFC)),
    n_above_05  = sum(abs_logFC > 0.5),
    n_total     = dplyr::n(),
    .groups     = "drop"
  ) |>
  dplyr::mutate(ci_lo = vapply(ci, `[[`, numeric(1), "lower"),
                ci_hi = vapply(ci, `[[`, numeric(1), "upper")) |>
  dplyr::select(-ci)

# ── Blunting ratio ----

med_tr_hr <- lfc_stats$med_abs_lfc[lfc_stats$contrast == "Training_HR"]
med_tr_lr <- lfc_stats$med_abs_lfc[lfc_stats$contrast == "Training_LR"]
blunting_ratio <- med_tr_hr / max(med_tr_lr, 1e-6)

# ── Cliff's delta pairwise ----

cliff_pairs <- list(
  c("Training_HR", "Training_LR"),
  c("Acute_HR",    "Acute_LR")
)

cliff_results <- purrr::map_dfr(cliff_pairs, function(pair) {
  x <- lfc_long$abs_logFC[lfc_long$contrast == pair[1]]
  y <- lfc_long$abs_logFC[lfc_long$contrast == pair[2]]
  tibble::tibble(contrast_1 = pair[1], contrast_2 = pair[2],
                 cliff_d = cliffs_delta(x, y))
})

# ── KS tests ----

ks_results <- purrr::map_dfr(cliff_pairs, function(pair) {
  x <- lfc_long$abs_logFC[lfc_long$contrast == pair[1]]
  y <- lfc_long$abs_logFC[lfc_long$contrast == pair[2]]
  ks <- ks.test(x, y)
  tibble::tibble(contrast_1 = pair[1], contrast_2 = pair[2],
                 ks_D = ks$statistic, ks_p = ks$p.value)
})

# ── Helper: build one histogram panel ----

st <- scale_text(BASE_STAT, PB_W)

make_hist <- function(ctr, show_x = TRUE) {
  d <- lfc_long |> dplyr::filter(contrast == ctr)
  s <- lfc_stats |> dplyr::filter(contrast == ctr)
  med_val <- s$med_lfc

  # Stat annotation (top-right, small)
  stat_label <- sprintf("med|FC| = %.2f [%.2f, %.2f]\nn(|FC|>0.5) = %d",
                        s$med_abs_lfc, s$ci_lo, s$ci_hi,
                        s$n_above_05)

  p <- ggplot(d, aes(x = logFC)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60,
                   fill = CONTRAST_COLORS[ctr],
                   color = "grey50", linewidth = 0.1, alpha = 0.80) +
    geom_density(linewidth = 0.5, color = darken_color(CONTRAST_COLORS[ctr], 0.5)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40",
               linewidth = 0.3) +
    geom_vline(xintercept = med_val, linetype = "dashed",
               color = darken_color(CONTRAST_COLORS[ctr], 0.5),
               linewidth = 0.6) +
    # Contrast title (top-left)
    annotate("text", x = -0.98, y = Inf, label = CTR_FACET[ctr],
             hjust = 0, vjust = 1.3, fontface = "bold",
             size = st * 1.05, color = "grey10") +
    # Stats below title (left-aligned, small italic)
    annotate("text", x = -0.98, y = Inf, label = stat_label,
             hjust = 0, vjust = 2.8, fontface = "italic",
             size = st * 0.68, color = "grey45", lineheight = 0.85) +
    # Median value label (right of dashed line, near top)
    annotate("label", x = med_val + ifelse(med_val >= 0, 0.04, -0.04),
             y = Inf,
             label = sprintf("%.3f", med_val),
             hjust = ifelse(med_val >= 0, 0, 1), vjust = 1.6,
             size = st * 0.65, fontface = "bold",
             color = darken_color(CONTRAST_COLORS[ctr], 0.5),
             fill = alpha("white", 0.85), label.size = 0,
             label.padding = unit(0.8, "pt")) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5),
                       expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    FIG_THEME +
    theme(legend.position  = "none",
          axis.text.y      = element_blank(),
          axis.ticks.y     = element_blank(),
          axis.title       = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin      = margin(1, 3, 1, 3))

  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }
  p
}

# ── Build individual panels ----

p_base  <- make_hist("Baseline_HRvLR", show_x = FALSE)
p_tr_hr <- make_hist("Training_HR",    show_x = FALSE)
p_tr_lr <- make_hist("Training_LR",    show_x = FALSE)
p_ac_hr <- make_hist("Acute_HR",       show_x = TRUE)
p_ac_lr <- make_hist("Acute_LR",       show_x = TRUE)

# ── Compose uniform 3×2 grid with spacer ----

# Condensed single-line subtitle
cliff_vals <- purrr::map_chr(seq_len(nrow(cliff_results)), function(i) {
  r <- cliff_results[i, ]
  tag <- if (grepl("Training", r$contrast_1)) "Tr" else "Ac"
  sprintf("%s=%.2f", tag, r$cliff_d)
})
ks_vals <- purrr::map_chr(seq_len(nrow(ks_results)), function(i) {
  r <- ks_results[i, ]
  tag <- if (grepl("Training", r$contrast_1)) "Tr" else "Ac"
  p_str <- if (r$ks_p < 0.001) "p<.001" else sprintf("p=%.3f", r$ks_p)
  sprintf("%s D=%.2f %s", tag, r$ks_D, p_str)
})
stat_sub <- sprintf("Blunt.=%.2f | Cliff d: %s | KS: %s",
                    blunting_ratio,
                    paste(cliff_vals, collapse = ", "),
                    paste(ks_vals, collapse = ", "))

pB_inner <- (p_base + p_tr_hr +
             p_tr_lr + p_ac_hr +
             p_ac_lr + plot_spacer()) +
  plot_layout(ncol = 2, nrow = 3) +
  plot_annotation(
    title    = expression(bold("Effect Size Distribution  ") *
                          "(dashed = median " * log[2] * "FC)"),
    subtitle = stat_sub,
    theme = theme(
      plot.title    = element_text(face = "bold", size = 9, hjust = 0),
      plot.subtitle = element_text(face = "italic", size = 6.5,
                                   color = "grey40", hjust = 0)
    )
  )

pB <- wrap_elements(pB_inner) +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE))

# ── Export audit tables ----

readr::write_csv(lfc_stats,
                 file.path(DAT_DIR, "audit_panel_B_median_lfc_ci.csv"))
readr::write_csv(cliff_results,
                 file.path(DAT_DIR, "audit_panel_B_cliff_delta.csv"))
readr::write_csv(ks_results,
                 file.path(DAT_DIR, "audit_panel_B_ks_tests.csv"))

# ── Save panel ----

save_panel(pB, file.path(RPT_DIR, "panel_B_logfc_density"),
           width = PB_W, height = PB_H)
