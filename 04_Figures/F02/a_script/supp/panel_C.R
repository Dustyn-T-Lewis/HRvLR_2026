# F02 Supp Panel C: responder heterogeneity (multivariate dispersion) + per-protein CV
# Reframes replicate variability from QC to biology: do HR and LR differ in how spread
# out their proteomes are? PERMDISP (Anderson 2006) tests homogeneity of multivariate
# dispersion per timepoint (each subject's distance to its group centroid); the CV
# scatter (Brenes 2019, linear scale) is the per-protein detail.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, patchwork, vegan)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

supp_dir <- file.path(RPT_DIR, "supp")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

PC_W <- 220
PC_H <- PANEL_H * 1.7

timepoints <- c("T1", "T2", "T3")

disp <- lapply(timepoints, function(tp) {
  ids <- meta$Col_ID[meta$Timepoint == tp & meta$Col_ID %in% colnames(imp_mat)]
  grp <- droplevels(meta$Group[match(ids, meta$Col_ID)])
  bd <- betadisper(dist(t(imp_mat[, ids])), grp)
  set.seed(42)
  p <- permutest(bd, permutations = 999)$tab$`Pr(>F)`[1]
  list(
    dist = tibble(timepoint = tp, Group = grp, dist_centroid = bd$distances),
    p = tibble(timepoint = tp, label = sprintf("PERMDISP %s", fmt_p(p)))
  )
})
disp_df <- bind_rows(lapply(disp, `[[`, "dist")) |>
  mutate(timepoint = factor(timepoint, levels = timepoints))
disp_p <- bind_rows(lapply(disp, `[[`, "p")) |>
  mutate(timepoint = factor(timepoint, levels = timepoints))

pC_disp <- ggplot(disp_df, aes(Group, dist_centroid, color = Group, fill = Group)) +
  geom_boxplot(alpha = 0.18, outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.12, height = 0, size = 1.4, alpha = 0.8) +
  geom_text(
    data = disp_p, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
    vjust = 1.4, size = FIG_GEOM_TEXT, fontface = "bold", color = "grey25"
  ) +
  facet_wrap(~timepoint) +
  scale_color_manual(values = GROUP_COLORS) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(
    title = "Responder Heterogeneity (multivariate dispersion)",
    subtitle = "distance to group centroid per subject | PERMDISP per timepoint",
    x = NULL, y = "distance to centroid"
  ) +
  FIG_THEME +
  theme(legend.position = "none")

lin_mat <- 2^as.matrix(imp_mat)
group_cv <- function(ids) {
  ids <- intersect(ids, colnames(lin_mat))
  if (length(ids) < 2) {
    return(rep(NA_real_, nrow(lin_mat)))
  }
  apply(lin_mat[, ids, drop = FALSE], 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) NA_real_ else sd(x) / mean(x) * 100
  })
}

cv_df <- lapply(timepoints, function(tp) {
  tibble(
    uniprot_id = imp_df$uniprot_id, gene = imp_df$gene, timepoint = tp,
    cv_HR = group_cv(meta$Col_ID[meta$Group == "HR" & meta$Timepoint == tp]),
    cv_LR = group_cv(meta$Col_ID[meta$Group == "LR" & meta$Timepoint == tp])
  )
}) |>
  bind_rows() |>
  filter(!is.na(cv_HR), !is.na(cv_LR)) |>
  mutate(timepoint = factor(timepoint, levels = timepoints), max_cv = pmax(cv_HR, cv_LR))

axis_cap <- as.numeric(quantile(c(cv_df$cv_HR, cv_df$cv_LR), 0.99, na.rm = TRUE))
cv_df$max_cv_capped <- pmin(cv_df$max_cv, as.numeric(quantile(cv_df$max_cv, 0.98, na.rm = TRUE)))

r_stats <- cv_df |>
  group_by(timepoint) |>
  summarise(n = n(), r = cor(cv_HR, cv_LR, use = "complete.obs"), .groups = "drop") |>
  rowwise() |>
  mutate(label = sprintf(
    "r = %.2f [%.2f, %.2f]", r,
    fisher_z_ci(r, n)[["lo"]], fisher_z_ci(r, n)[["hi"]]
  )) |>
  ungroup()

pC_cv <- ggplot(cv_df, aes(cv_HR, cv_LR)) +
  facet_wrap(~timepoint) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = max_cv_capped), alpha = 0.35, size = 0.7) +
  geom_label(
    data = r_stats, aes(x = 0, y = axis_cap, label = label), inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = FIG_GEOM_TEXT, color = "grey30", fontface = "bold",
    fill = scales::alpha("white", 0.85), label.size = 0, label.padding = unit(1.5, "pt")
  ) +
  scale_color_viridis_c(
    option = "inferno", direction = -1, begin = 0.1, end = 0.85, name = "CV%",
    guide = guide_colorbar(barwidth = unit(2, "mm"), barheight = unit(12, "mm"))
  ) +
  coord_equal(xlim = c(0, axis_cap), ylim = c(0, axis_cap)) +
  labs(
    title = "Per-Protein Variability (CV%)",
    subtitle = "linear-scale CV per protein, HR vs LR within timepoint (Brenes 2019)",
    x = expression(bold(CV * "%"[HR])), y = expression(bold(CV * "%"[LR]))
  ) +
  FIG_THEME +
  theme(legend.position = "right", legend.key.size = unit(3, "mm"))

pC <- (pC_disp / pC_cv) +
  plot_annotation(
    tag_levels = list(c("C", "")),
    theme = theme(plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE))
  )

save_png(pC, file.path(supp_dir, "panel_c_variability"), PC_W, PC_H)
write.csv(disp_df, file.path(DAT_DIR, "audit_panel_C_dispersion.csv"), row.names = FALSE)
write.csv(cv_df |> select(uniprot_id, gene, timepoint, cv_HR, cv_LR),
  file.path(DAT_DIR, "audit_panel_C_cv_scatter.csv"),
  row.names = FALSE
)
cat("F02 Supp Panel C done.\n")
