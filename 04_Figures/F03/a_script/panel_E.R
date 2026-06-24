# F03 Panel E: linkage verdict (permutation null + sPLS CV) ----
setwd(rprojroot::find_rstudio_root_file())
if (!exists("perm_null")) source("04_Figures/F03/a_script/HRvLR_F03_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Pull Q2 values from spls_cv ----
q2_df <- spls_cv |>
  dplyr::filter(metric == "Q2") |>
  dplyr::select(component, q2 = value) |>
  dplyr::mutate(component = factor(component, levels = sort(unique(component))))

obs_max_t <- perm_null$observed_max_t[1]
p_fam <- perm_null$p_perm[1]

# Left: permutation null ----
p_null <- ggplot(perm_null, aes(x = null_max_t)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(xintercept = obs_max_t, colour = "#B2182B", linewidth = 1) +
  annotate(
    "text",
    x = obs_max_t, y = Inf,
    label = sprintf(
      "observed max|t| = %.2f\nfamily-wise p = %.3f",
      obs_max_t, p_fam
    ),
    hjust = 1.05, vjust = 1.4, size = 3, colour = "#B2182B"
  ) +
  labs(
    x = "Max |t| under subject-permuted phenotype",
    y = "Permutations",
    title = NULL
  ) +
  FIG_THEME

# Right: sPLS cross-validated Q2 ----
q2_min <- min(q2_df$q2)
y_lo <- min(q2_min * 1.08, q2_min - 0.3)

p_spls <- ggplot(q2_df, aes(x = component, y = q2)) +
  geom_col(fill = "#4393C3", colour = "black", linewidth = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate(
    "text",
    x = 1.5, y = y_lo * 0.55,
    label = "Q2 <= 0: does not generalise",
    size = 3, colour = "grey30", hjust = 0.5
  ) +
  scale_y_continuous(limits = c(y_lo, max(q2_df$q2) + 0.05)) +
  labs(x = "sPLS component", y = "Cross-validated Q2") +
  FIG_THEME

# Combine ----
panel_e <- (p_null | p_spls) +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(
    title = "The module-phenotype linkage does not survive validation",
    subtitle = "Left: strongest grid signal vs permuted null (p=0.096). Right: sPLS Q2<=0 out-of-sample. Null at n~16.",
    tag_levels = list(c("E", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE, color = "grey30"),
      plot.tag      = element_text(face = "bold", size = FIG_TAG_SIZE)
    )
  )

save_panel(panel_e, file.path(RPT_DIR, "panel_e_verdict"), 200, 90)
write.csv(q2_df, file.path(DAT_DIR, "audit_panel_E_spls_q2.csv"), row.names = FALSE)
cat("F03 Panel E done.\n")
