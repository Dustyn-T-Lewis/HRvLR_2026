# 04_Figures - HRvLR Unified Style
# Single source of truth: palettes, themes, sizing constants, helpers.

pacman::p_load(ggplot2, scales, grid)

# Colour palettes. Two blue/red mappings coexist and must not be conflated:
# GROUP_COLORS encode responder (HR = dark blue, LR = dark red); DIR_COLORS encode
# logFC direction (Up = light red, Down = light blue). The group hues are the
# darker shades so a responder legend never reads as a direction legend. These are
# HRvLR-specific: the YvO/Mito suites use the lighter #4393C3/#D6604D for their own
# groups, so blue/red meaning does not carry across suites - always read the legend.
GROUP_COLORS <- c(HR = "#2166AC", LR = "#B2182B")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_FILL <- c(
  HR_T1 = scales::alpha("#2166AC", 0.40),
  HR_T2 = "#2166AC",
  HR_T3 = scales::alpha("#2166AC", 0.70),
  LR_T1 = scales::alpha("#B2182B", 0.40),
  LR_T2 = "#B2182B",
  LR_T3 = scales::alpha("#B2182B", 0.70)
)

# Timepoint hues (Okabe-Ito, colourblind-safe): baseline / trained / acute
TIME_COLORS <- c(T1 = "#E69F00", T2 = "#0072B2", T3 = "#009E73")

# Contrast codes: responder family sets the hue (HR blue / LR red / between-responder
# green / interaction purple); timepoint sets the shade (training lighter -> acute darker,
# baseline lightest). One palette for every contrast-keyed panel.
CONTRAST_COLORS <- c(
  Training_HR          = "#6BAED6", # HR, training (T2-T1)
  Acute_HR             = "#2166AC", # HR, acute (T3-T2)
  Training_LR          = "#FC9272", # LR, training
  Acute_LR             = "#B2182B", # LR, acute
  Baseline_HRvLR       = "#A1D99B", # HR vs LR, baseline (T1)
  Trained_HRvLR        = "#41AB5D", # HR vs LR, trained (T2)
  Acute_HRvLR          = "#238B45", # HR vs LR, acute (T3)
  Training_Interaction = "#9E9AC8", # differential training response
  Acute_Interaction    = "#6A51A3" # differential acute response
)

# Phenotype domains for F01: muscle size, fibre size, strength. Reuses the
# established hues so the palette keeps a single source.
DOMAIN_COLORS <- c(
  fibre = unname(CONTRAST_COLORS["Training_HR"]),
  muscle = unname(GROUP_COLORS["HR"]),
  strength = "grey50"
)

# GO ontology hues for enrichment panels (BP / CC / MF).
ONT_COLORS <- c(BP = "#66C2A5", CC = "#FC8D62", MF = "#8DA0CB")

# Group x timepoint shades for the 6-cell ordination: HR blues, LR reds, each
# ramping light (T1) to dark (T3) within its hue.
GROUP_TIME_COLORS <- c(
  HR_T1 = "#9ECAE1", HR_T2 = "#4292C6", HR_T3 = "#08519C",
  LR_T1 = "#FCAE91", LR_T2 = "#EF3B2C", LR_T3 = "#99000D"
)

# Uniform render height (mm) for every F02 panel so they tile at one scale.
PANEL_H <- 95

# Shared in-plot text size (mm) for counts, stat annotations, and point labels.
FIG_GEOM_TEXT <- 2.6

# Theme

FIG_BASE_SIZE <- 10
FIG_TITLE_SIZE <- 12
FIG_SUBTITLE_SIZE <- 9
FIG_STRIP_SIZE <- 10
FIG_AXIS_TITLE <- 10
FIG_AXIS_TEXT <- 8.5
FIG_LEGEND_TITLE <- 9.5
FIG_LEGEND_TEXT <- 8.5
FIG_TAG_SIZE <- 15

FIG_THEME <- theme_bw(base_size = FIG_BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle = element_text(
      face = "bold.italic", size = FIG_SUBTITLE_SIZE,
      color = "grey30"
    ),
    plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title = element_text(face = "bold", size = FIG_AXIS_TITLE),
    axis.text = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title = element_text(
      face = "bold", size = FIG_LEGEND_TITLE,
      color = "grey20"
    ),
    legend.text = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size = unit(3, "mm"),
    panel.grid.minor = element_blank()
  )

# Key/legend sizing

KEY_TEXT <- 2.8

# Contrast label mappings

CTR_SHORT <- c(
  Training_HR          = "Tr.(HR)",
  Training_LR          = "Tr.(LR)",
  Acute_HR             = "Ac.(HR)",
  Acute_LR             = "Ac.(LR)",
  Baseline_HRvLR       = "Base.",
  Trained_HRvLR        = "Trn.(HRvLR)",
  Acute_HRvLR          = "Ac.(HRvLR)",
  Training_Interaction = "Tr.Int.",
  Acute_Interaction    = "Ac.Int."
)

# Helper functions

fmt_p <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  }
  if (p < 0.01) {
    return(sprintf("p = %.3f", p))
  }
  sprintf("p = %.2f", p)
}

fisher_z_ci <- function(r, n, k = 0, level = 0.95) {
  n_eff <- n - k
  if (n_eff < 4 || is.na(r)) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  z <- atanh(r)
  se <- 1 / sqrt(n_eff - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * se), hi = tanh(z + crit * se))
}

# Bonett & Wright 2000, Psychometrika 65(1), doi:10.1007/bf02294183: the improved
# standard error for a Spearman correlation. Use this, not fisher_z_ci, whenever r
# is a Spearman rho.
fisher_z_ci_spearman <- function(r, n, level = 0.95) {
  if (n < 4 || is.na(r)) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  z <- atanh(r)
  se <- sqrt((1 + r^2 / 2) / (n - 3))
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * se), hi = tanh(z + crit * se))
}

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ ""
  )
}

clean_pathway_name <- function(name, max_chars = 45) {
  name |>
    stringr::str_remove("^HALLMARK_") |>
    stringr::str_remove("^GOBP_") |>
    stringr::str_remove("^GOCC_") |>
    stringr::str_remove("^GOMF_") |>
    stringr::str_remove("^REACTOME_") |>
    stringr::str_remove("^KEGG_MEDICUS_") |>
    stringr::str_remove("^KEGG_") |>
    stringr::str_remove("^GOSLIM_") |>
    stringr::str_remove("^WP_") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_title() |>
    stringr::str_replace("Mtorc1", "mTORC1") |>
    stringr::str_replace("Myc ", "MYC ") |>
    stringr::str_replace("E2f ", "E2F ") |>
    stringr::str_replace("Dna ", "DNA ") |>
    stringr::str_replace("Rna ", "RNA ") |>
    stringr::str_replace("Tnfa ", "TNFa ") |>
    stringr::str_replace("Uv ", "UV ") |>
    stringr::str_replace("G2m ", "G2M ") |>
    stringr::str_replace("Il6 ", "IL6 ") |>
    stringr::str_replace("Il2 ", "IL2 ") |>
    stringr::str_replace("Kras ", "KRAS ") |>
    stringr::str_replace("P53 ", "p53 ") |>
    stringr::str_replace("Tgf ", "TGF ") |>
    stringr::str_replace("Nf Kb", "NF-kB") |>
    stringr::str_replace("Atp ", "ATP ") |>
    stringr::str_replace("Nadh ", "NADH ") |>
    stringr::str_replace("Oxidative Phosphorylation", "OXPHOS") |>
    stringr::str_trunc(max_chars, ellipsis = "...")
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun, ...)
}

scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

# Floating bracket above tallest data point (pad = fraction of range)
bracket_pos <- function(y, pad = 0.08) {
  y <- y[!is.na(y)]
  max(y) + pad * diff(range(y))
}

# Device & export helpers

get_pdf_device <- function() {
  tryCatch(
    {
      cairo_pdf(tempfile())
      dev.off()
      cairo_pdf
    },
    error = function(e) "pdf"
  )
}

# Cache once at source() time - avoids 12+ temp-file probes
PDF_DEVICE <- get_pdf_device()

save_panel <- function(plot, path_stem, width, height, pdf_device = PDF_DEVICE) {
  ggsave(paste0(path_stem, ".pdf"), plot,
    width = width, height = height, units = "mm", device = pdf_device
  )
  ggsave(paste0(path_stem, ".png"), plot,
    width = width, height = height, units = "mm", dpi = 300
  )
}

save_png <- function(plot, path_stem, width, height) {
  ggsave(paste0(path_stem, ".png"), plot,
    width = width, height = height, units = "mm", dpi = 300
  )
}

# F03 database palette

DB_COLORS <- c(
  Hallmark = "#AA336A", KEGG = "#E65100",
  Reactome = "#1565C0", "GO:BP" = "#00796B",
  "GO:CC" = "#0097A7", "GO:MF" = "#7CB342",
  "GO Slim" = "#5D4037"
)

# F03 ring-volcano palettes, one per contrast family so a reader never conflates
# a within-group trajectory with a between-responder difference or an interaction.
# nes is the diverging arc ramp (down pole, midpoint, up pole); up/down recolour
# the volcano points to the same poles. Responses read as the canonical red/blue;
# baseline is green/purple; the interaction takes orange and its own purple.
RING_PALETTES <- list(
  responses = list(nes = c("#2166AC", "grey95", "#B2182B"), up = "#B2182B", down = "#2166AC"),
  differential = list(nes = c("#6A51A3", "grey95", "#238B45"), up = "#238B45", down = "#6A51A3"),
  interaction = list(nes = c("#762A83", "grey95", "#E65100"), up = "#E65100", down = "#762A83")
)
