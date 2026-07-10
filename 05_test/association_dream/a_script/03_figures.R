# Panels for the dream association fits: an effect-size volcano of the features
# whose trajectory tracks the continuous phenotype, a timepoint heatmap of the top
# phenotype-associated pathways, and the categorical HR-vs-LR effect summary.
if (!exists("scores")) source(here::here("05_test", "association_dream", "a_script", "setup.R"))
pacman::p_load(stringr, ggrepel, forcats)

cont_res <- read_csv(file.path(dat, "dream_continuous.csv"), show_col_types = FALSE)
cat_res <- read_csv(file.path(dat, "dream_categorical.csv"), show_col_types = FALSE)

trait_lab <- function(x) unname(cont_traits[x])
coef_cols <- c(
  phenotype = "#7B1FA2", `phenotype:T2` = "#0072B2", `phenotype:T3` = "#009E73"
)

volc_label <- function(d, n = 8) {
  d |>
    group_by(matrix) |>
    slice_min(p_value, n = n, with_ties = FALSE) |>
    ungroup()
}

fig_volcano <- ggplot(
  cont_res |> mutate(matrix = recode(matrix, protein = "Proteins", pathway = "Pathways")),
  aes(estimate, -log10(p_value), color = coefficient)
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(size = 1.1, alpha = 0.55) +
  ggrepel::geom_text_repel(
    seed = 42,
    data = volc_label(cont_res) |>
      mutate(
        matrix = recode(matrix, protein = "Proteins", pathway = "Pathways"),
        lab = paste0(
          ifelse(matrix == "Pathways", clean_pathway_name(feature, 22), feature),
          " (", recode(phenotype,
            comp_hypertrophy = "comp", d_fcsa_I = "fCSA-I", d_fcsa_II = "fCSA-II",
            d_mcsa = "mCSA", d_1rm_legpress = "1RM-leg", d_1rm_ext = "1RM-ext"
          ), ")"
        )
      ),
    aes(label = lab), size = 2.1, color = "grey20",
    max.overlaps = 20, min.segment.length = 0, box.padding = 0.4
  ) +
  scale_color_manual(values = coef_cols, name = "coefficient") +
  facet_wrap(~matrix, scales = "free_x") +
  labs(
    x = "dream coefficient (per SD of phenotype)", y = "-log10 p",
    title = "Features whose trajectory tracks the continuous phenotype",
    subtitle = "dream: feature ~ phenotype * timepoint + ACCUM_VL + (1|subject); no feature survives BH FDR 0.05"
  ) +
  FIG_THEME
save_panel(fig_volcano, file.path(rpt, "fig1_continuous_volcano"), 220, 130)

top_pw <- cont_res |>
  filter(matrix == "pathway") |>
  slice_min(p_value, n = 15, with_ties = FALSE) |>
  distinct(feature) |>
  pull(feature)

hm <- scores[top_pw, ] |>
  as.data.frame() |>
  rownames_to_column("pathway") |>
  pivot_longer(-pathway, names_to = "sample", values_to = "score") |>
  left_join(meta, by = "sample") |>
  group_by(pathway, timepoint) |>
  summarise(mean_score = mean(score), .groups = "drop") |>
  group_by(pathway) |>
  mutate(z = as.numeric(scale(mean_score))) |>
  ungroup() |>
  mutate(pathway = fct_reorder(clean_pathway_name(pathway, 34), z, .fun = function(x) x[1]))

fig_heat <- ggplot(hm, aes(timepoint, pathway, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = unname(DIR_COLORS["Down"]), mid = "grey95",
    high = unname(DIR_COLORS["Up"]), midpoint = 0,
    name = "row z"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Top phenotype-associated pathways across timepoints",
    subtitle = "Mean singscore per timepoint, z-scored within pathway; baseline T1 to acute T3"
  ) +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 7), panel.grid = element_blank())
save_panel(fig_heat, file.path(rpt, "fig2_pathway_heatmap"), 170, 150)

cat_top <- cat_res |>
  mutate(matrix = recode(matrix, protein = "Proteins", pathway = "Pathways")) |>
  group_by(matrix, coefficient) |>
  slice_min(p_value, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    lab = ifelse(matrix == "Pathways", clean_pathway_name(feature, 28), feature),
    term = reorder_within(lab, -log10(p_value), interaction(matrix, coefficient)),
    dir = ifelse(estimate > 0, "higher in HR", "higher in LR")
  )

fig_cat <- ggplot(cat_top, aes(-log10(p_value), term, color = dir)) +
  geom_segment(aes(x = 0, xend = -log10(p_value), yend = term), linewidth = 0.5) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "grey55") +
  scale_color_manual(values = c(
    "higher in HR" = unname(GROUP_COLORS["HR"]),
    "higher in LR" = unname(GROUP_COLORS["LR"])
  ), name = NULL) +
  scale_y_reordered() +
  facet_wrap(matrix ~ coefficient, scales = "free_y", ncol = 3) +
  labs(
    x = "-log10 p", y = NULL,
    title = "HR-vs-LR categorical effects (group and group x timepoint)",
    subtitle = "dream: feature ~ group * timepoint + ACCUM_VL + (1|subject); dotted = nominal p = 0.05"
  ) +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 6))
save_panel(fig_cat, file.path(rpt, "fig3_categorical_effects"), 230, 170)

message("association_dream figures written")
