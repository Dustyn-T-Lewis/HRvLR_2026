#!/usr/bin/env Rscript
# The mCSA axis in four panels: what the phenotype is made of, how its parts are
# weighted, what a scan of the 931 complete-case proteins against d_mcsa
# returns, and what the one survivor looks like.
#
# Each panel title states a claim rather than naming its contents. No caption
# block sits on the render; the caveats live in F_mcsa_axis_legend.md beside it.
#
# Read the sheet left to right, top to bottom: the two CSA axes separate but one
# subject sets how far; the composite's formula and the fitted factor disagree
# about how much whole-muscle CSA matters; only the T2 scan leaves its null; and
# the protein that does it is concurrent, not a baseline forecast.

pacman::p_load(
  here, dplyr, tidyr, readr, tibble, ggplot2, ggrepel, patchwork, forcats
)
source(here("functions", "shared_style.R"))
source(here("functions", "shared_hlm.R"))
source(here("functions", "blood_index_model.R"))
source(here("03_Features", "04_galamm_pilot", "a_script", "pilot_helpers.R"))
source(here("03_Features", "06_mcsa_axis", "a_script", "mcsa_helpers.R"))

DAT <- here("03_Features", "06_mcsa_axis", "c_data")
RPT <- here("03_Features", "06_mcsa_axis", "b_reports")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

SET_LINES <- c(`all subjects` = "solid", `LR_S14 dropped` = "22")
CONFIG_LABEL <- c(
  T1 = "T1 baseline", T2 = "T2 trained", T3 = "T3 acute",
  delta = "T2 - T1 delta"
)
CSA_LABEL <- c(
  d_fcsa_I = "fCSA I", d_fcsa_II = "fCSA II", d_mcsa = "whole-muscle CSA"
)

pheno <- phenotype_table()
read_cell <- function(f) read_csv(file.path(DAT, f), show_col_types = FALSE)
weights <- read_cell("01_composite_weights.csv")
scans <- read_cell("02_protein_mcsa.csv")
survivors <- read_cell("02_survivors.csv")

symbols <- read_csv(
  here("03_Features", "01_Proteins", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
) |>
  select(feature = "uniprot_id", gene = "gene")

axes <- pheno |>
  mutate(fibre = fibre_axis(pheno), arm = .data$group_arm)
trimmed <- filter(axes, .data$subject != DISCORDANT_SUBJECT)
rho_all <- stats::cor(axes$d_mcsa, axes$fibre, method = "spearman")
rho_cut <- stats::cor(trimmed$d_mcsa, trimmed$fibre, method = "spearman")

p_axes <- ggplot(axes, aes(.data$d_mcsa, .data$fibre, colour = .data$arm)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 2.4) +
  geom_text_repel(
    data = filter(axes, .data$subject == DISCORDANT_SUBJECT),
    aes(label = .data$subject), size = FIG_GEOM_TEXT, seed = 42,
    nudge_x = -1, nudge_y = -0.35, show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf, y = -Inf, hjust = -0.08, vjust = -0.5, size = FIG_GEOM_TEXT,
    colour = "grey25", lineheight = 1.1,
    label = sprintf(
      "rho = %.2f all 16\nrho = %.2f without LR_S14", rho_all, rho_cut
    )
  ) +
  scale_colour_manual(values = GROUP_COLORS, name = NULL) +
  labs(
    title = "Fibre and whole-muscle growth are two axes",
    subtitle = "and one subject sets how far apart they sit",
    x = "Whole-muscle CSA change (d_mcsa)",
    y = "Fibre axis (mean z of fCSA I, II)", tag = "A"
  ) +
  FIG_THEME +
  theme(legend.position = c(0.13, 0.9), legend.background = element_blank())

shares <- weights |>
  transmute(
    item = .data$item,
    `Composite formula` = .data$weight / sum(.data$weight),
    `Fitted factor` = .data$loading / sum(.data$loading)
  ) |>
  pivot_longer(-"item", names_to = "weighting", values_to = "share") |>
  mutate(item = fct_rev(factor(CSA_LABEL[.data$item], levels = CSA_LABEL)))

p_weights <- ggplot(
  shares, aes(.data$share, .data$item, colour = .data$weighting)
) +
  geom_line(aes(group = .data$item), colour = "grey75", linewidth = 0.6) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(
    `Composite formula` = "grey35", `Fitted factor` = DIR_COLORS[["Up"]]
  ), name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.55)) +
  labs(
    title = "The two weightings disagree about mCSA",
    subtitle = "share of the three CSA measures under each",
    x = "Share of the weighting", y = NULL, tag = "B"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

qq <- scans |>
  mutate(config = factor(CONFIG_LABEL[.data$config], levels = CONFIG_LABEL)) |>
  arrange(.data$config, .data$subject_set, .data$p) |>
  mutate(
    expected = -log10(seq_along(.data$p) / (n() + 1)),
    observed = -log10(.data$p),
    .by = c("config", "subject_set")
  )
bh_line <- scans |>
  filter(.data$bh < 0.05) |>
  summarise(cut = -log10(max(.data$p)))
top_hit <- qq |>
  semi_join(
    mutate(survivors, config = CONFIG_LABEL[.data$config]),
    by = c("feature", "subject_set", "config")
  ) |>
  left_join(symbols, by = "feature")

p_scan <- ggplot(
  qq, aes(.data$expected, .data$observed, colour = .data$config)
) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(
    yintercept = bh_line$cut, linetype = "dotted",
    colour = DIR_COLORS[["Up"]], linewidth = 0.4
  ) +
  geom_line(aes(linetype = .data$subject_set), linewidth = 0.55) +
  geom_point(data = top_hit, size = 1.8, show.legend = FALSE) +
  geom_text_repel(
    data = top_hit, aes(label = .data$gene), size = FIG_GEOM_TEXT,
    seed = 42, nudge_x = -0.5, nudge_y = 0.4, show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 0, y = bh_line$cut, hjust = -0.05, vjust = -0.6,
    label = "BH q = 0.05", size = FIG_GEOM_TEXT - 0.2,
    colour = DIR_COLORS[["Up"]]
  ) +
  scale_colour_manual(
    values = unname(c(TIME_COLORS, delta = "grey40")), name = NULL
  ) +
  scale_linetype_manual(values = SET_LINES, name = NULL) +
  labs(
    title = "Only the trained timepoint leaves its null",
    subtitle = "931 proteins per config; only T1 precedes the outcome",
    x = expression(Expected ~ -log[10] ~ italic(p)),
    y = expression(Observed ~ -log[10] ~ italic(p)), tag = "C"
  ) +
  FIG_THEME +
  theme(legend.position = "right", legend.spacing.y = unit(1, "mm"))

hit <- survivors$feature[1]
hit_gene <- symbols$gene[match(hit, symbols$feature)]
inp <- pilot_data()
hit_long <- c("T1", "T2", "T3") |>
  lapply(function(tp) {
    design <- config_design(inp$mat, inp$meta, pheno, tp)
    tibble(
      config = tp, subject = design$subject, d_mcsa = design$mcsa,
      abundance = design$y[hit, ]
    )
  }) |>
  bind_rows() |>
  mutate(
    arm = ifelse(startsWith(.data$subject, "HR"), "HR", "LR"),
    config = factor(CONFIG_LABEL[.data$config], levels = CONFIG_LABEL)
  )
hit_q <- scans |>
  filter(
    .data$feature == hit, .data$subject_set == "all subjects",
    .data$config != "delta"
  ) |>
  transmute(
    config = factor(CONFIG_LABEL[.data$config], levels = CONFIG_LABEL),
    label = sprintf("q = %.3f", .data$bh)
  )

p_hit <- ggplot(hit_long, aes(.data$d_mcsa, .data$abundance)) +
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE, colour = "grey30",
    fill = "grey85", linewidth = 0.5
  ) +
  geom_point(aes(colour = .data$arm), size = 2.1) +
  geom_text(
    data = hit_q, aes(x = -Inf, y = Inf, label = .data$label),
    hjust = -0.15, vjust = 1.4, size = FIG_GEOM_TEXT, fontface = "bold",
    colour = "grey20", inherit.aes = FALSE
  ) +
  facet_wrap(~config) +
  scale_colour_manual(values = GROUP_COLORS, name = NULL) +
  labs(
    title = sprintf("%s tracks growth after training, not before", hit_gene),
    subtitle = "same protein, same subjects; T1 is the only forecast",
    x = "Whole-muscle CSA change (d_mcsa)",
    y = sprintf("%s abundance (log2)", hit_gene), tag = "D"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

sheet <- (p_axes | p_weights) / (p_scan | p_hit)
save_panel(sheet, file.path(RPT, "F_mcsa_axis"), 260, 200)
cat("F_mcsa_axis written\n")
