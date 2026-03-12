# Figure 1 — Panel D: DEPs per Contrast (Pseudo-log Stacked Bar) ----
# Exports: sig_sets, dir_map, all_genes, pi_total, fdr_total (used by Panel E)

if (!exists("meta")) source("04_Figures/F1/a_script/HRvLR_F1_setup.R")

RPT <- RPT_DIR
DAT <- DAT_DIR

PD_W <- 170

SET_LABELS <- CTR_FACET[DISPLAY_ORDER]
all_genes <- unique(dep_df$gene[!is.na(dep_df$gene)])

sig_sets <- list()
dir_map  <- list()

for (ctr in CONTRASTS) {
  pi_vals  <- dep_df[[paste0("pi_score_", ctr)]]
  lfc_vals <- dep_df[[paste0("logFC_", ctr)]]
  is_sig   <- !is.na(pi_vals) & pi_vals < 0.05
  sig_sets[[ctr]] <- dep_df$gene[is_sig]
  dir_map[[ctr]]  <- setNames(ifelse(lfc_vals[is_sig] > 0, "Up", "Down"),
                               dep_df$gene[is_sig])
}

n_total <- length(all_genes)
pi_ci <- data.frame(
  contrast = CONTRASTS,
  n_sig    = vapply(sig_sets, length, integer(1)),
  n_total  = n_total,
  pct      = 100 * vapply(sig_sets, length, integer(1)) / n_total,
  ci_lo    = vapply(sig_sets, function(s) 100 * binom.test(length(s), n_total)$conf.int[1], double(1)),
  ci_hi    = vapply(sig_sets, function(s) 100 * binom.test(length(s), n_total)$conf.int[2], double(1))
)

pi_total  <- sum(vapply(sig_sets, length, integer(1)))
fdr_total <- sum(vapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  if (fdr_col %in% names(dep_df)) sum(dep_df[[fdr_col]] < 0.05, na.rm = TRUE) else 0
}, integer(1)))

frac_list <- lapply(DISPLAY_ORDER, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  p_col   <- paste0("P.Value_", ctr)
  tibble(
    contrast  = SET_LABELS[ctr],
    threshold = c("p < 0.05", "q < 0.05", "\u03A0 < 0.05"),
    n = c(sum(!is.na(dep_df[[p_col]])   & dep_df[[p_col]]   < 0.05),
          sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.05),
          length(sig_sets[[ctr]]))
  )
})
frac_df <- bind_rows(frac_list) |>
  mutate(
    contrast  = factor(contrast, levels = rev(unname(SET_LABELS))),
    threshold = factor(threshold, levels = c("p < 0.05", "q < 0.05", "\u03A0 < 0.05")),
    pct       = 100 * n / length(all_genes),
    fill_key  = paste(contrast, threshold, sep = "___")
  ) |>
  filter(n > 1)

SET_DISPLAY_COLORS <- setNames(unname(CONTRAST_COLORS[DISPLAY_ORDER]),
                                unname(SET_LABELS))

FRAC_FILL <- c()
for (cname in names(SET_DISPLAY_COLORS)) {
  col <- unname(SET_DISPLAY_COLORS[cname])
  FRAC_FILL[paste(cname, "p < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.25)
  FRAC_FILL[paste(cname, "q < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.55)
  FRAC_FILL[paste(cname, "\u03A0 < 0.05", sep = "___")] <- col
}

THRESH_LABEL <- c("p < 0.05" = "p \u2264 0.05", "q < 0.05" = "FDR \u2264 0.05",
                  "\u03A0 < 0.05" = "\u03A0 \u2264 0.05")

label_df <- frac_df |>
  group_by(contrast) |> arrange(contrast, threshold) |>
  mutate(label     = THRESH_LABEL[as.character(threshold)],
         next_pct  = lead(pct, default = 0),
         seg_width = pct - next_pct,
         label_y   = (next_pct + pct) / 2,
         text_col  = if_else(threshold == "p < 0.05", "grey20", "white")) |>
  filter(seg_width > 0.3) |>
  ungroup()

# Plot ----

n_ctr <- length(DISPLAY_ORDER)
bg_rects <- lapply(seq_len(n_ctr), function(i) {
  annotate("rect",
           xmin = n_ctr - i + 0.5, xmax = n_ctr - i + 1.5,
           ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS[DISPLAY_ORDER[i]], alpha = 0.20,
           color = "grey70", linewidth = 0.2)
})

pD <- ggplot(frac_df, aes(x = contrast, y = pct, fill = fill_key))
for (r in bg_rects) pD <- pD + r
pD <- pD +
  geom_col(position = "identity", width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(data = label_df,
            aes(x = contrast, y = label_y, label = label, color = I(text_col)),
            inherit.aes = FALSE, hjust = 0.5,
            size = scale_text(BASE_COUNT - 1.0, PD_W), fontface = "bold") +
  scale_fill_manual(values = FRAC_FILL) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = exp(1)),
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  coord_flip() +
  labs(title = "DEPs per Contrast",
       subtitle = sprintf("Fraction of %s filtered proteins | \u03A0: %d / FDR: %d total",
                          format(length(all_genes), big.mark = ","),
                          pi_total, fdr_total),
       x = NULL, y = "% of proteome (pseudo-log scale)",
       tag = "D") +
  FIG_THEME + theme(legend.position = "none",
                    axis.text.y = element_text(face = "bold", size = 6.5))

# Export ----

write_csv(pi_ci, file.path(DAT, "audit_panel_D_dep_fraction_ci.csv"))

save_panel(pD, file.path(RPT, "panel_D_dep_counts"),
           width = PD_W, height = 90)

# sig_sets, dir_map, all_genes, pi_total, fdr_total are now available
# in the calling environment (composite sources this script at top level)
