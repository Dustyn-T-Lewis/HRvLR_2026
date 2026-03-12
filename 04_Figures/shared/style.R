# 04_Figures — HRvLR Unified Style ----
# Single source of truth: palettes, themes, sizing constants, helpers.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(grid)
})

# ── Colour palettes ----

GROUP_COLORS <- c(HR = "#2166AC", LR = "#B2182B")
DIR_COLORS   <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_FILL <- c(
  HR_T1 = scales::alpha("#2166AC", 0.40),
  HR_T2 = "#2166AC",
  HR_T3 = scales::alpha("#2166AC", 0.70),
  LR_T1 = scales::alpha("#B2182B", 0.40),
  LR_T2 = "#B2182B",
  LR_T3 = scales::alpha("#B2182B", 0.70)
)

SHAPE_TP <- c(T1 = 16, T2 = 17, T3 = 15)

CONTRAST_COLORS <- c(
  Training_HR          = "#2166AC",
  Training_LR          = "#B2182B",
  Acute_HR             = "#67A9CF",
  Acute_LR             = "#EF8A62",
  Baseline_HRvLR       = "#4CAF50",
  Training_Interaction = "#9B7FBF",
  Acute_Interaction    = "#FF8F00"
)

PCA_COLORS <- c(
  HR_T1 = scales::alpha("#2166AC", 0.45),
  HR_T2 = "#2166AC",
  HR_T3 = "#4393C3",
  LR_T1 = scales::alpha("#B2182B", 0.45),
  LR_T2 = "#B2182B",
  LR_T3 = "#D6604D"
)

PCA_SHAPES <- c(HR_T1 = 16, HR_T2 = 17, HR_T3 = 15,
                LR_T1 = 16, LR_T2 = 17, LR_T3 = 15)

# ── Panel dimensions & sizing ----

PANEL_SM <- 120
PANEL_MD <- 180
PANEL_LG <- 280

BASE_PATHWAY  <- 4.0
BASE_GENE     <- 3.2
BASE_STAT     <- 3.5
BASE_QUADRANT <- 4.0
BASE_COUNT    <- 3.5
BASE_TAG      <- 18

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * sqrt(panel_width_mm / ref_width)
}

# ── Theme ----

FIG_BASE_SIZE     <- 10
FIG_TITLE_SIZE    <- 12
FIG_SUBTITLE_SIZE <- 9
FIG_STRIP_SIZE    <- 10
FIG_AXIS_TITLE    <- 10
FIG_AXIS_TEXT     <- 8.5
FIG_LEGEND_TITLE  <- 9.5
FIG_LEGEND_TEXT   <- 8.5
FIG_TAG_SIZE      <- 15

FIG_THEME <- theme_bw(base_size = FIG_BASE_SIZE) +
  theme(
    plot.title       = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle    = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                    color = "grey30"),
    plot.tag         = element_text(face = "bold", size = FIG_TAG_SIZE),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title       = element_text(face = "bold", size = FIG_AXIS_TITLE),
    axis.text        = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title     = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                    color = "grey20"),
    legend.text      = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size  = unit(3, "mm"),
    panel.grid.minor = element_blank()
  )

# ── Key/legend sizing ----

KEY_TEXT <- 2.8

# ── Contrast definitions ----

CONTRASTS <- c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR",
               "Baseline_HRvLR", "Training_Interaction", "Acute_Interaction")

# ── Contrast label mappings ----

CTR_FACET <- c(
  Training_HR          = "Training (HR)",
  Training_LR          = "Training (LR)",
  Acute_HR             = "Acute (HR)",
  Acute_LR             = "Acute (LR)",
  Baseline_HRvLR       = "Baseline (HR vs LR)",
  Training_Interaction = "Training Interaction",
  Acute_Interaction    = "Acute Interaction"
)

CTR_SHORT <- c(
  Training_HR          = "Tr.(HR)",
  Training_LR          = "Tr.(LR)",
  Acute_HR             = "Ac.(HR)",
  Acute_LR             = "Ac.(LR)",
  Baseline_HRvLR       = "Base.",
  Training_Interaction = "Tr.Int.",
  Acute_Interaction    = "Ac.Int."
)

# ── Helper functions ----

fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01)  return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

fisher_z_ci <- function(r, n, k = 0, level = 0.95) {
  n_eff <- n - k
  if (n_eff < 4 || is.na(r)) return(c(lo = NA_real_, hi = NA_real_))
  z  <- atanh(r)
  se <- 1 / sqrt(n_eff - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * se), hi = tanh(z + crit * se))
}

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ ""
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

boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

# Floating bracket above tallest data point (pad = fraction of range)
bracket_pos <- function(y, pad = 0.08) {
  y <- y[!is.na(y)]
  max(y) + pad * diff(range(y))
}

darken_color <- function(col, factor = 0.7) {
  rgb_vals <- grDevices::col2rgb(col) / 255
  vapply(seq_along(col), function(i)
    grDevices::rgb(rgb_vals[1, i] * factor, rgb_vals[2, i] * factor,
                   rgb_vals[3, i] * factor),
    character(1))
}

strip_plot_meta <- function(p) {
  p + theme(plot.title = element_blank(), plot.subtitle = element_blank())
}

# ── Device & export helpers ----

get_pdf_device <- function() {
  tryCatch(
    { cairo_pdf(tempfile()); dev.off(); cairo_pdf },
    error = function(e) "pdf"
  )
}

# Cache once at source() time — avoids 12+ temp-file probes
PDF_DEVICE <- get_pdf_device()

save_panel <- function(plot, path_stem, width, height, pdf_device = PDF_DEVICE) {
  ggsave(paste0(path_stem, ".pdf"), plot,
         width = width, height = height, units = "mm", device = pdf_device)
  ggsave(paste0(path_stem, ".png"), plot,
         width = width, height = height, units = "mm", dpi = 300)
}
