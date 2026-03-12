# Figure 1 — Panel A: CV% Violins (Inter-Individual Variability) ----

if (!exists("meta")) source("04_Figures/F1/a_script/HRvLR_F1_setup.R")

PA_W <- 160
PA_H <- 100

# Computation ----

lin_mat <- 2^as.matrix(norm_df[, samp_names])

cv_list <- lapply(levels(meta$group), function(g) {
  idx <- meta$sample_id[meta$group == g]
  sub <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
  tibble(group = g, cv = cv_pct)
})
cv_df <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group <- factor(cv_df$group,
                      levels = c("HR_T1", "HR_T2", "HR_T3",
                                 "LR_T1", "LR_T2", "LR_T3"))

# Bootstrap 95% CI on median CV per group
set.seed(42)

ci_df <- cv_df |>
  group_by(group) |>
  summarise(
    med = median(cv),
    ci  = list(boot_median_ci(cv)),
    .groups = "drop"
  ) |>
  mutate(ci_lo = vapply(ci, `[[`, numeric(1), "lower"),
         ci_hi = vapply(ci, `[[`, numeric(1), "upper")) |>
  select(-ci)

cv_med <- ci_df |> select(group, med)

# Statistical tests: Wilcoxon + Cliff's delta + BH correction ----

bracket_comps <- list(
  c("HR_T1", "LR_T1"),   # baseline HR vs LR
  c("HR_T1", "HR_T2"),   # training HR
  c("LR_T1", "LR_T2"),   # training LR
  c("HR_T1", "HR_T3"),   # acute HR
  c("LR_T1", "LR_T3")    # acute LR
)
bracket_pvals_raw <- vapply(bracket_comps, function(pair)
  wilcox.test(cv_df$cv[cv_df$group == pair[1]],
              cv_df$cv[cv_df$group == pair[2]])$p.value, double(1))
bracket_pvals <- p.adjust(bracket_pvals_raw, method = "BH")

cliff_results <- data.frame(
  comparison = vapply(bracket_comps, paste, character(1), collapse = " vs "),
  p_raw      = bracket_pvals_raw,
  p_bh       = bracket_pvals,
  cliffs_d   = vapply(bracket_comps, function(pair)
    cliffs_delta(cv_df$cv[cv_df$group == pair[1]],
                 cv_df$cv[cv_df$group == pair[2]]), double(1))
)

# Filter for significance AND non-negligible effect (|d| >= 0.15)
sig_idx <- which(bracket_pvals < 0.05 & abs(cliff_results$cliffs_d) >= 0.15)

cv_ymax      <- quantile(cv_df$cv, 0.99)
bracket_base <- cv_ymax * 0.85
bracket_step <- cv_ymax * 0.08

# Plot ----

group_labels <- c(HR_T1 = "HR T1", HR_T2 = "HR T2", HR_T3 = "HR T3",
                  LR_T1 = "LR T1", LR_T2 = "LR T2", LR_T3 = "LR T3")

pA <- ggplot(cv_df, aes(x = group, y = cv, fill = group)) +
  geom_violin(alpha = 0.7, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, linewidth = 0.3, fill = "white") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_label(data = cv_med, aes(x = group, y = med, label = sprintf("%.0f%%", med)),
             size = scale_text(BASE_COUNT - 0.5, PA_W),
             fontface = "bold", fill = alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1, "pt")) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_x_discrete(labels = group_labels) +
  labs(title = "Inter-Individual Variability (CV%)",
       subtitle = "Protein-level CV on cycloess-normalized intensities (linear scale)",
       x = NULL, y = "CV (%)", tag = "A") +
  FIG_THEME + theme(legend.position = "none")

.pd_label <- function(p, d) {
  p_str <- if (p < 0.001) "p < .001" else sprintf("p = %.3f", p)
  sprintf("%s\nd = %.2f", p_str, abs(d))
}

pd_sz   <- scale_text(KEY_TEXT - 1.2, PA_W)
brk_col <- "grey25"

if (length(sig_idx) > 0) {
  sig_ypos <- bracket_base + (seq_along(sig_idx) - 1) * bracket_step
  tip <- cv_ymax * 0.015
  for (i in seq_along(sig_idx)) {
    pair <- bracket_comps[[sig_idx[i]]]
    x1 <- match(pair[1], levels(cv_df$group))
    x2 <- match(pair[2], levels(cv_df$group))
    yb <- sig_ypos[i]
    xmid <- (x1 + x2) / 2
    pA <- pA +
      annotate("segment", x = x1, xend = x2, y = yb, yend = yb,
               linewidth = 0.5, color = brk_col) +
      annotate("segment", x = x1, xend = x1, y = yb - tip, yend = yb,
               linewidth = 0.5, color = brk_col) +
      annotate("segment", x = x2, xend = x2, y = yb - tip, yend = yb,
               linewidth = 0.5, color = brk_col) +
      annotate("label", x = xmid, y = yb + cv_ymax * 0.01,
               label = .pd_label(bracket_pvals[sig_idx[i]],
                                  cliff_results$cliffs_d[sig_idx[i]]),
               size = pd_sz, fontface = "bold", color = brk_col,
               fill = alpha("white", 0.85), linewidth = 0.15,
               label.padding = unit(1.2, "pt"),
               hjust = 0.5, vjust = 0.5, lineheight = 0.85)
  }
  y_upper <- max(sig_ypos) + bracket_step * 2.5
  pA <- pA + coord_cartesian(ylim = c(0, y_upper), clip = "on")
} else {
  pA <- pA + coord_cartesian(ylim = c(0, bracket_base + bracket_step),
                              clip = "on")
}

# Export ----

write_csv(ci_df, file.path(DAT_DIR, "audit_panel_A_median_cv_ci.csv"))
write_csv(cliff_results, file.path(DAT_DIR, "audit_panel_A_wilcoxon_effects.csv"))

save_panel(pA, file.path(RPT_DIR, "panel_A_cv"),
           width = PA_W, height = PA_H)
