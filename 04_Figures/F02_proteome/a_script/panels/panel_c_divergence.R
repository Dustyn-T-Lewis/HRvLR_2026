# F02 Panel C: responder divergence leads, protein level (HR vs LR logFC)
# Per protein, the training/acute logFC in HR against the same in LR. On-diagonal =
# concordant response; off-diagonal = divergent. Candidates come from the interaction
# contrast pi-gate (sig_pi on Training/Acute_Interaction), the right form of the test
# (Gelman & Stern 2006). The pi-gate is not FDR-controlled and none survive BH, so
# these are exploratory leads; the panel prints both counts.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, ggrepel)

if (!exists("meta")) source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))

PE_W <- 120
PE_H <- 110

phase_levels <- c("Training (T2 - T1)", "Acute (T3 - T2)")

diverg_df <- lapply(c(Training = "Training", Acute = "Acute"), function(ph) {
  tibble(
    gene = dep_df$gene,
    hr = dep_df[[paste0("logFC_", ph, "_HR")]],
    lr = dep_df[[paste0("logFC_", ph, "_LR")]],
    sig_int = dep_df[[paste0("sig_pi_", ph, "_Interaction")]],
    fdr = dep_df[[paste0("adj.P.Val_", ph, "_Interaction")]]
  ) |>
    filter(!is.na(gene), !is.na(hr), !is.na(lr)) |>
    mutate(
      phase = ph,
      divergent = !is.na(sig_int) & sig_int != 0,
      fdr_sig = !is.na(fdr) & fdr < 0.05,
      gap = abs(hr - lr)
    )
}) |>
  bind_rows() |>
  mutate(phase = factor(recode(phase, Training = phase_levels[1], Acute = phase_levels[2]), levels = phase_levels))

# Label the top 5 up- and top 5 down-interaction leads per phase (direction from
# the HR-minus-LR sign), so both arms of the divergence are named.
label_df <- diverg_df |>
  filter(divergent) |>
  mutate(direction = if_else(hr - lr > 0, "up", "down")) |>
  group_by(phase, direction) |>
  slice_max(gap, n = 5, with_ties = FALSE) |>
  ungroup()

count_df <- diverg_df |>
  group_by(phase) |>
  summarise(n_div = sum(divergent), n_fdr = sum(fdr_sig), .groups = "drop") |>
  mutate(label = sprintf("%d pi-gated / %d at FDR", n_div, n_fdr))

# Each facet gets its own tight symmetric range so one outlier can't stretch both.
blank_df <- diverg_df |>
  filter(divergent) |>
  group_by(phase) |>
  summarise(lim = max(abs(c(hr, lr))) * 1.18, .groups = "drop") |>
  tidyr::crossing(sgn = c(-1, 1)) |>
  transmute(phase, hr = sgn * lim, lr = sgn * lim)
INT_COLOR <- unname(CONTRAST_COLORS[["Acute_Interaction"]])

pC <- ggplot(diverg_df, aes(hr, lr)) +
  facet_wrap(~phase, scales = "free") +
  geom_blank(data = blank_df) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_point(data = ~ filter(.x, !divergent), color = "grey75", alpha = 0.3, size = 0.7) +
  geom_point(data = ~ filter(.x, divergent), color = INT_COLOR, alpha = 0.9, size = 1.4) +
  geom_text_repel(
    data = label_df, aes(label = gene), size = FIG_GEOM_TEXT, fontface = "bold",
    color = INT_COLOR, min.segment.length = 0, segment.size = 0.2,
    box.padding = 0.6, point.padding = 0.25, max.overlaps = Inf, seed = 42
  ) +
  geom_label(
    data = count_df, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -0.06, vjust = 1.4, size = FIG_GEOM_TEXT, fontface = "bold", color = INT_COLOR,
    fill = scales::alpha("white", 0.85), label.size = 0, label.padding = unit(1.5, "pt")
  ) +
  labs(
    x = expression(bold(log[2] * FC[HR])), y = expression(bold(log[2] * FC[LR]))
  ) +
  FIG_THEME

save_png(pC, file.path(RPT_DIR, "panels", "panel_c_divergence"), PE_W, PE_H)
F02_AUDIT[["panel_C_divergence"]] <- diverg_df |> filter(divergent) |> select(phase, gene, hr, lr, gap)
cat("F02 Panel C done.\n")

# --- Supplement (owned by this panel): responder heterogeneity.
# PERMDISP on multivariate dispersion + the per-protein CV scatter behind it.
pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, patchwork, vegan)


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
F02_AUDIT[["panel_C_dispersion"]] <- disp_df
F02_AUDIT[["panel_C_cv_scatter"]] <- cv_df |> select(uniprot_id, gene, timepoint, cv_HR, cv_LR)
cat("F02 Supp Panel C done.\n")
