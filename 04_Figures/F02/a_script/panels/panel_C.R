# F02 Panel C: Per-protein CV% HR vs LR scatter, one facet per timepoint
# Inter-individual CV% (linear scale, Brenes 2024) compared between responder
# groups within each timepoint. Replaces the redundant CV violin / transition pair.

pacman::p_load(here, dplyr, tidyr, ggplot2)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

PC_W <- 200
PC_H <- PANEL_H

# CV on the linear scale per protein within a set of samples
lin_mat <- 2^as.matrix(imp_mat)

group_cv <- function(ids) {
  ids <- intersect(ids, colnames(lin_mat))
  if (length(ids) < 2) {
    return(rep(NA_real_, nrow(lin_mat)))
  }
  apply(lin_mat[, ids, drop = FALSE], 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) {
      return(NA_real_)
    }
    sd(x) / mean(x) * 100
  })
}

timepoints <- c("T1", "T2", "T3")
cv_df <- lapply(timepoints, function(tp) {
  hr_ids <- meta$Col_ID[meta$Group == "HR" & meta$Timepoint == tp]
  lr_ids <- meta$Col_ID[meta$Group == "LR" & meta$Timepoint == tp]
  tibble(
    uniprot_id = imp_df$uniprot_id,
    gene = imp_df$gene,
    timepoint = tp,
    cv_HR = group_cv(hr_ids),
    cv_LR = group_cv(lr_ids)
  )
}) |>
  bind_rows() |>
  filter(!is.na(cv_HR), !is.na(cv_LR)) |>
  mutate(
    timepoint = factor(timepoint, levels = timepoints),
    max_cv = pmax(cv_HR, cv_LR)
  )

cv_cap <- as.numeric(quantile(cv_df$max_cv, 0.98, na.rm = TRUE))
cv_df$max_cv_capped <- pmin(cv_df$max_cv, cv_cap)
axis_cap <- as.numeric(quantile(c(cv_df$cv_HR, cv_df$cv_LR), 0.99, na.rm = TRUE))

r_stats <- cv_df |>
  group_by(timepoint) |>
  summarise(
    n = sum(!is.na(cv_HR) & !is.na(cv_LR)),
    r = cor(cv_HR, cv_LR, use = "complete.obs"),
    .groups = "drop"
  ) |>
  rowwise() |>
  mutate(
    ci_lo = fisher_z_ci(r, n)[["lo"]],
    ci_hi = fisher_z_ci(r, n)[["hi"]],
    label = sprintf("r = %.2f [%.2f, %.2f]", r, ci_lo, ci_hi)
  ) |>
  ungroup()

pC <- ggplot(cv_df, aes(cv_HR, cv_LR)) +
  facet_wrap(~timepoint, nrow = 1) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "grey50", linewidth = 0.4
  ) +
  geom_point(aes(color = max_cv_capped), alpha = 0.35, size = 0.7) +
  geom_label(
    data = r_stats, aes(x = 0, y = axis_cap, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = FIG_GEOM_TEXT,
    color = "grey30", fontface = "bold", fill = scales::alpha("white", 0.85),
    label.size = 0, label.padding = unit(1.5, "pt")
  ) +
  scale_color_viridis_c(
    option = "inferno", direction = -1, begin = 0.1, end = 0.85,
    name = "CV%",
    guide = guide_colorbar(barwidth = unit(2, "mm"), barheight = unit(12, "mm"))
  ) +
  coord_equal(xlim = c(0, axis_cap), ylim = c(0, axis_cap)) +
  labs(
    title = "Per-Protein Variability (CV%)",
    subtitle = sprintf(
      "linear-scale CV per protein, HR vs LR within timepoint (Brenes 2024) | r: T1 %.2f, T2 %.2f, T3 %.2f",
      r_stats$r[r_stats$timepoint == "T1"],
      r_stats$r[r_stats$timepoint == "T2"],
      r_stats$r[r_stats$timepoint == "T3"]
    ),
    x = expression(bold(CV * "%"[HR])),
    y = expression(bold(CV * "%"[LR])),
    tag = "C"
  ) +
  FIG_THEME +
  theme(
    legend.position = "right",
    legend.key.size = unit(3, "mm")
  )

save_png(pC, file.path(RPT_DIR, "panels", "panel_c_cv_scatter"), PC_W, PC_H)
write.csv(cv_df |> select(uniprot_id, gene, timepoint, cv_HR, cv_LR),
  file.path(DAT_DIR, "audit_panel_C_cv_scatter.csv"),
  row.names = FALSE
)
cat("F02 Panel C done.\n")
