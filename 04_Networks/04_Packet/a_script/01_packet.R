# The network packet: how the modules were built, what each one is, how
# they move, then the contrasts and the two screens.
# Reads only 04_Networks c_data. Pages describe; none of them concludes.

pacman::p_load(here, dplyr, tidyr, forcats, purrr, ggplot2, patchwork)

source(here("functions", "contrasts.R"))
source(here("functions", "screen_pages.R"))

stage <- function(...) here("04_Networks", ...)
modules <- readRDS(stage("01_Modules", "c_data", "modules.rds"))
char <- readRDS(stage(
  "02_Characterise", "c_data", "module_characterisation.rds"
))
screens <- readRDS(stage(
  "03_Classify_Associate", "c_data", "module_screens.rds"
))
module_colours <- rownames(modules$eigengenes) |> (\(m) setNames(m, m))()
module_levels <- char$labels |>
  arrange(desc(.data$n_proteins)) |>
  pull("module")

sft <- as_tibble(modules$soft_threshold) |>
  mutate(signed_r2 = -sign(.data$slope) * .data$SFT.R.sq)
sft_panel <- function(y, label) {
  ggplot(sft, aes(.data$Power, .data[[y]])) +
    geom_line(colour = "grey60") +
    geom_point(aes(colour = .data$Power == modules$power), size = 2) +
    scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#B2182B")) +
    labs(x = "Soft power", y = label) +
    FIG_THEME +
    theme(legend.position = "none")
}
p_sft <- sft_panel("signed_r2", "Signed scale-free R2") +
  geom_hline(yintercept = 0.85, linetype = "dashed") +
  sft_panel("mean.k.", "Mean connectivity") +
  plot_annotation(
    title = "Soft-threshold choice",
    subtitle = sprintf(
      "WGCNA signed, bicor, subject-centred; power %d (R2 %.2f, k %.1f)",
      modules$power, modules$r2, modules$mean_k
    ),
    caption = caption(
      "Left: scale-free fit at each power; dashed line is pickSoftThreshold's ",
      "0.85 criterion, and the red point is the lowest power that clears it. ",
      "Right: mean connectivity at each power. ",
      "Data: 04_Networks/01_Modules/c_data/modules.xlsx, sheet soft_threshold."
    ),
    theme = FIG_THEME
  )

size_icc <- char$labels |>
  left_join(modules$icc, by = "module") |>
  mutate(module = factor(.data$module, levels = rev(module_levels)))
p_sizes <- ggplot(size_icc, aes(.data$n_proteins, .data$module)) +
  geom_col(aes(fill = .data$module), colour = "grey30", linewidth = 0.2) +
  geom_text(
    aes(label = sprintf("ICC %.2f", .data$icc)),
    hjust = -0.15, size = 3
  ) +
  scale_fill_manual(values = module_colours, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title = "Module sizes and subject dependence",
    subtitle = sprintf(
      "%d modules; %d of %d proteins unassigned (grey)",
      nrow(modules$eigengenes), sum(modules$membership$module == "grey"),
      nrow(modules$membership)
    ),
    x = "Proteins", y = NULL,
    caption = caption(
      "Bars: proteins per module. Label: intraclass correlation of the ",
      "eigengene across a subject's three biopsies (lme4, 1 | subject). ",
      "Modules are defined on within-subject-centred data so that subject ",
      "identity cannot build them; a high ICC here means the scored ",
      "eigengene still differs between people. ",
      "Hubs (highest kME): ",
      paste(
        sprintf("%s: %s", char$labels$module, char$labels$hubs),
        collapse = "; "
      ),
      ". Data: 04_Networks/01_Modules/c_data/modules.xlsx, sheet subject_icc."
    )
  ) +
  FIG_THEME

p_ora <- char$ora_fit |>
  clusterProfiler::filter(.data$p.adjust < 0.05) |>
  enrichplot::dotplot(showCategory = 3, label_format = 50, font.size = 7) +
  labs(
    title = "What each module is: enrichment",
    subtitle = sprintf(
      "clusterProfiler ORA, universe %d; top 3 sets per module, BH < 0.05",
      nrow(modules$membership)
    ),
    caption = caption(
      "enrichplot::dotplot of the compareCluster result. Columns: modules ",
      "(members with a set in brackets). Dot: one of the three sets with the ",
      "smallest p in that module at BH < 0.05; size = gene ratio, colour = ",
      "BH within module; modules with none are absent. Hubs per module: ",
      "caption of the module-size page. ",
      "Data: 04_Networks/02_Characterise/c_data/module_characterisation.xlsx, ",
      "sheet ora."
    )
  ) +
  theme(
    plot.caption = element_text(hjust = 0, size = 8, colour = "grey25"),
    plot.caption.position = "plot"
  )

p_string <- char$string |>
  mutate(module = factor(.data$module, levels = rev(module_levels))) |>
  ggplot(aes(y = .data$module)) +
  geom_segment(
    aes(x = .data$expected, xend = .data$edges, yend = .data$module),
    colour = "grey70"
  ) +
  geom_point(aes(x = .data$expected), shape = 124, size = 4) +
  geom_point(aes(x = .data$edges, fill = .data$module), shape = 21, size = 3) +
  geom_text(
    aes(x = .data$edges, label = sprintf("%.1fx", .data$ratio)),
    hjust = -0.4, size = 2.8
  ) +
  scale_fill_manual(values = module_colours, guide = "none") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "STRING edges within each module",
    subtitle = sprintf(
      "STRINGdb PPI enrichment, v12, score >= %d; detected background",
      char$string_min_score
    ),
    x = "Within-module edges (log scale)", y = NULL,
    caption = caption(
      "Point: observed high-confidence STRING edges among a module's ",
      "members. Tick: edges STRINGdb expects from the members' degrees in ",
      "the detected background. Label: observed / expected. Modules at ",
      "BH < 0.05: ", sum(char$string$bh < 0.05), " of ", nrow(char$string),
      ". Data: 04_Networks/02_Characterise/c_data/",
      "module_characterisation.xlsx, sheet string."
    )
  ) +
  FIG_THEME

me_long <- as_tibble(t(modules$eigengenes), rownames = "sample_id") |>
  pivot_longer(-"sample_id", names_to = "module", values_to = "me") |>
  left_join(modules$meta, by = "sample_id") |>
  mutate(module = factor(.data$module, levels = module_levels))
p_traj <- me_long |>
  summarise(
    mean = mean(.data$me), se = stats::sd(.data$me) / sqrt(n()),
    .by = c("module", "arm", "timepoint")
  ) |>
  ggplot(aes(
    .data$timepoint, .data$mean,
    colour = .data$arm, group = .data$arm
  )) +
  geom_line(
    data = me_long, aes(y = .data$me, group = .data$subject),
    alpha = 0.2, linewidth = 0.3
  ) +
  geom_line(linewidth = 0.8) +
  geom_pointrange(
    aes(ymin = .data$mean - .data$se, ymax = .data$mean + .data$se),
    size = 0.25
  ) +
  facet_wrap(~module, nrow = 3) +
  scale_colour_manual(values = GROUP_COLORS) +
  labs(
    title = "Eigengene trajectories",
    subtitle = "Module eigengene per biopsy, scored on raw abundance",
    x = NULL, y = "Eigengene", colour = "Arm",
    caption = caption(
      "Thin lines: one subject's three biopsies. Thick line and range: arm ",
      "mean +/- standard error at each timepoint. ",
      "Data: 04_Networks/01_Modules/c_data/modules.xlsx, sheet eigengenes."
    )
  ) +
  FIG_THEME

p_contrasts <- screens$contrasts |>
  mutate(
    contrast = factor(.data$contrast, levels = CONTRAST_NAMES),
    module = factor(.data$module, levels = rev(module_levels))
  ) |>
  ggplot(aes(.data$contrast, .data$module, fill = .data$t)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(.data$p < 0.05, sprintf("%.3f", .data$p), "")),
    size = 2.6
  ) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  labs(
    title = "Module eigengenes across the nine contrasts",
    subtitle = sprintf(
      "limma, protein design, subject block (cor %.3f); BH per contrast",
      screens$correlation$consensus
    ),
    x = NULL, y = NULL, fill = "Moderated t",
    caption = caption(
      "Fill: moderated t of each module in each contrast. Label: nominal p ",
      "where below 0.05. Eigengenes at BH < 0.05 across all contrasts: ",
      sum(screens$contrasts$bh < 0.05), ". ",
      "Data: 04_Networks/03_Classify_Associate/c_data/module_screens.xlsx, ",
      "sheet contrasts."
    )
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

SCREENS <- "04_Networks/03_Classify_Associate/c_data/module_screens.xlsx"
p_auc <- auc_page(screens$classify, "Module", SCREENS)
p_assoc <- association_page(screens$chance_associate, "Module", SCREENS)

write_packet(
  c(
    list(
      "Soft-threshold choice" = p_sft,
      "Module sizes and subject dependence" = p_sizes,
      "What each module is: enrichment" = p_ora,
      "STRING edges within each module" = p_string,
      "Eigengene trajectories" = p_traj,
      "Module eigengenes across the nine contrasts" = p_contrasts,
      "Module AUC per classification task" = p_auc,
      "Module association with phenotype" = p_assoc
    ),
    screen_hit_pages(screens, "Module", SCREENS)
  ),
  stage("04_Packet", "b_reports", "04_Networks_packet.pdf"),
  title = "HRvLR 04 Networks"
)
