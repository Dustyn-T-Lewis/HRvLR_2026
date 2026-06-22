# F02 Panel D: CV% Comparison (Pre vs Post Normalization) ----
# Faceted by Group (HR | LR), Pre/Post x-axis. Median labels with bootstrap CIs.
# Adapted from YvO F02 panel_D pattern.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

PD_W <- 110
PD_H <- 120
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

# All six Group_Time cells across the 2x3 design
GROUP_ORDER <- c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
GROUP_COLS <- GROUP_FILL

# CV on linear (not log) scale per group
lin_mat <- 2^as.matrix(imp_mat)

cv_list <- lapply(GROUP_ORDER, function(g) {
  idx <- meta$Col_ID[meta$Group_Time == g]
  idx <- intersect(idx, colnames(lin_mat))
  if (length(idx) < 2) {
    return(NULL)
  }
  sub <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) {
      return(NA_real_)
    }
    sd(x) / mean(x) * 100
  })
  tibble(group = g, cv = cv_pct)
})
cv_df <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group <- factor(cv_df$group, levels = GROUP_ORDER)

cv_df$resp_group <- factor(ifelse(grepl("HR", cv_df$group), "HR", "LR"),
  levels = c("HR", "LR")
)
cv_df$timepoint <- factor(sub("^(HR|LR)_", "", as.character(cv_df$group)),
  levels = c("T1", "T2", "T3")
)

# Bootstrap 95% CI on median CV per group
set.seed(42)

# Pairwise Wilcoxon tests, BH corrected (audit only)
bracket_comps <- list(
  c("HR_T1", "HR_T2"), c("HR_T2", "HR_T3"),
  c("LR_T1", "LR_T2"), c("LR_T2", "LR_T3"),
  c("HR_T1", "LR_T1"), c("HR_T2", "LR_T2"), c("HR_T3", "LR_T3")
)
bracket_pvals_raw <- sapply(bracket_comps, function(pair) {
  wilcox.test(
    cv_df$cv[cv_df$group == pair[1]],
    cv_df$cv[cv_df$group == pair[2]]
  )$p.value
})
bracket_pvals <- p.adjust(bracket_pvals_raw, method = "BH")

cliff_results <- data.frame(
  comparison = sapply(bracket_comps, paste, collapse = " vs "),
  p_raw = bracket_pvals_raw,
  p_bh = bracket_pvals,
  cliffs_d = sapply(bracket_comps, function(pair) {
    cliffs_delta(
      cv_df$cv[cv_df$group == pair[1]],
      cv_df$cv[cv_df$group == pair[2]]
    )
  })
)

cv_ci <- cv_df |>
  group_by(resp_group, timepoint, group) |>
  summarise(
    med = median(cv),
    ci_lo = boot_median_ci(cv)[["lower"]],
    ci_hi = boot_median_ci(cv)[["upper"]],
    cv_max = max(cv),
    .groups = "drop"
  )

n_prot <- nrow(imp_df)
grand_med <- median(cv_df$cv)
grand_ci <- boot_median_ci(cv_df$cv)

# Median CV per timepoint, wide, for the two transition arrows
delta_cv <- cv_ci |>
  select(resp_group, timepoint, med) |>
  pivot_wider(names_from = timepoint, values_from = med)

# Two arrows per group: training (T1->T2) and acute (T2->T3)
arrow_df <- bind_rows(
  delta_cv |> transmute(resp_group,
    x = 1, xend = 2, y = T1, yend = T2,
    x_mid = 1.5, y_mid = (T1 + T2) / 2, arrow_label = sprintf("%+.1f%%", T2 - T1)
  ),
  delta_cv |> transmute(resp_group,
    x = 2, xend = 3, y = T2, yend = T3,
    x_mid = 2.5, y_mid = (T2 + T3) / 2, arrow_label = sprintf("%+.1f%%", T3 - T2)
  )
)

sub_txt <- sprintf(
  "%s proteins | CV %.0f%% [%.0f\u2013%.0f]\ntraining (T1\u2192T2) HR %+.1f%%, LR %+.1f%% | acute (T2\u2192T3) HR %+.1f%%, LR %+.1f%%",
  format(n_prot, big.mark = ","), grand_med, grand_ci[1], grand_ci[2],
  delta_cv$T2[delta_cv$resp_group == "HR"] - delta_cv$T1[delta_cv$resp_group == "HR"],
  delta_cv$T2[delta_cv$resp_group == "LR"] - delta_cv$T1[delta_cv$resp_group == "LR"],
  delta_cv$T3[delta_cv$resp_group == "HR"] - delta_cv$T2[delta_cv$resp_group == "HR"],
  delta_cv$T3[delta_cv$resp_group == "LR"] - delta_cv$T2[delta_cv$resp_group == "LR"]
)

pD <- ggplot(cv_df, aes(x = timepoint, y = cv, fill = group)) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_quasirandom(aes(color = group),
    alpha = 0.15, size = 0.5,
    width = 0.25, groupOnX = TRUE, show.legend = FALSE
  ) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA, linewidth = 0.3,
    color = "black", fill = "white", coef = 0
  ) +
  geom_hline(
    yintercept = 25, linetype = "dashed", color = "grey50",
    linewidth = 0.4
  ) +
  # Median labels with bootstrap CI
  geom_label(
    data = cv_ci,
    aes(
      x = timepoint, y = cv_max + 3,
      label = sprintf("%.0f%% [%.0f\u2013%.0f]", med, ci_lo, ci_hi)
    ),
    size = scale_text(BASE_COUNT + 0.5, PD_W / 2),
    fontface = "bold", fill = alpha("white", 0.8),
    linewidth = 0.2, label.padding = unit(1.5, "pt"),
    hjust = 0, nudge_x = -0.4
  ) +
  # Directional arrows: training (T1->T2) and acute (T2->T3) median shifts
  geom_segment(
    data = arrow_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, color = "grey30",
    arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
    linewidth = 0.6
  ) +
  geom_label(
    data = arrow_df,
    aes(x = x_mid, y = y_mid, label = arrow_label),
    inherit.aes = FALSE, size = scale_text(BASE_COUNT, PD_W / 2),
    fontface = "bold.italic", fill = alpha("white", 0.85),
    label.padding = unit(1.5, "pt"), linewidth = 0.2,
    color = "grey30"
  ) +
  facet_wrap(~resp_group, nrow = 1) +
  scale_fill_manual(values = GROUP_COLS) +
  scale_color_manual(values = GROUP_COLS) +
  coord_cartesian(ylim = c(0, max(cv_ci$cv_max) + 12)) +
  labs(
    title = "Inter-Individual Variability (CV%)",
    subtitle = sub_txt,
    x = NULL, y = "CV (%)",
    tag = "D"
  ) +
  FIG_THEME +
  theme(
    legend.position = "none",
    panel.spacing = unit(8, "mm")
  )

write.csv(as.data.frame(cv_ci),
  file.path(DAT_DIR, "audit_panel_D_median_cv_ci.csv"),
  row.names = FALSE
)
write.csv(cliff_results,
  file.path(DAT_DIR, "audit_panel_D_wilcoxon_effects.csv"),
  row.names = FALSE
)

save_panel(pD, file.path(RPT_DIR, "panel_d_cv_transitions"), PD_W, PD_H)
cat("F02 Panel D done.\n")
