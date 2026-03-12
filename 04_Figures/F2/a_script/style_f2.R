# style_f2.R — HRvLR F2 style definitions
# Adapts shared/style.R with F2-specific scatter classification and ring helpers.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
})

# ── Colour palettes ───────────────────────────────────────────────────────────

DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_COLORS <- c(HR = "#2166AC", LR = "#B2182B")

SIG_COLORS_SCATTER <- c(
  "HR only"     = "#E05A4E",
  "LR only"     = "#5DA5DA",
  "Sig Both"    = "#2E7D32",
  "Interaction" = "#7B5EA7",
  "NS"          = "grey70"
)

SIG_LABEL_FILL_SCATTER <- c(
  "HR only"     = scales::alpha("#E05A4E", 0.80),
  "LR only"     = scales::alpha("#5DA5DA", 0.80),
  "Sig Both"    = scales::alpha("#2E7D32", 0.80),
  "Interaction" = scales::alpha("#7B5EA7", 0.80),
  "NS"          = scales::alpha("grey70",  0.75)
)

SIG_LABEL_TEXT_SCATTER <- c(
  "HR only"     = "white",
  "LR only"     = "white",
  "Sig Both"    = "white",
  "Interaction" = "white",
  "NS"          = "white"
)

# ── Contrasts ─────────────────────────────────────────────────────────────────

CONTRASTS <- c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR",
               "Baseline_HRvLR", "Training_Interaction", "Acute_Interaction")

CTR_SHORT <- c(
  Training_HR           = "Trn_HR",
  Training_LR           = "Trn_LR",
  Acute_HR              = "Acu_HR",
  Acute_LR              = "Acu_LR",
  Baseline_HRvLR        = "Base",
  Training_Interaction  = "Trn_Int",
  Acute_Interaction     = "Acu_Int"
)

CONTRAST_COLORS <- c(
  Training_HR           = "#D6604D",
  Training_LR           = "#4393C3",
  Acute_HR              = "#B2182B",
  Acute_LR              = "#2166AC",
  Baseline_HRvLR        = "#74C476",
  Training_Interaction  = "#7B5EA7",
  Acute_Interaction     = "#E08214"
)

# ── Panel dimensions & sizing ─────────────────────────────────────────────────

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

# ── Theme ─────────────────────────────────────────────────────────────────────

FIG_THEME <- theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(face = "bold.italic", size = 9, color = "grey30"),
    plot.tag         = element_text(face = "bold", size = 15),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    axis.title       = element_text(face = "bold", size = 10),
    axis.text        = element_text(size = 8.5, color = "grey15"),
    legend.title     = element_text(face = "bold", size = 9.5, color = "grey20"),
    legend.text      = element_text(size = 8.5, color = "grey15"),
    legend.key.size  = unit(3, "mm"),
    panel.grid.minor = element_blank()
  )

# ── Scatter classification ────────────────────────────────────────────────────

classify_proteins_scatter <- function(pi_HR, pi_LR, pi_interaction,
                                       threshold = 0.05) {
  dplyr::case_when(
    pi_interaction < threshold             ~ "Interaction",
    pi_HR < threshold & pi_LR < threshold  ~ "Sig Both",
    pi_HR < threshold                      ~ "HR only",
    pi_LR < threshold                      ~ "LR only",
    TRUE                                   ~ "NS"
  ) |>
    factor(levels = c("Interaction", "Sig Both", "HR only", "LR only", "NS"))
}

# ── Helpers ───────────────────────────────────────────────────────────────────

get_pdf_device <- function() {
  tryCatch(
    { f <- tempfile(); cairo_pdf(f); dev.off(); unlink(f); cairo_pdf },
    error = function(e) "pdf"
  )
}

PDF_DEVICE <- get_pdf_device()

save_panel <- function(plot, path_stem, width, height, pdf_device = PDF_DEVICE) {
  ggsave(paste0(path_stem, ".pdf"), plot,
         width = width, height = height, units = "mm", device = pdf_device)
  ggsave(paste0(path_stem, ".png"), plot,
         width = width, height = height, units = "mm", dpi = 300)
  message("  Saved: ", basename(path_stem))
}
