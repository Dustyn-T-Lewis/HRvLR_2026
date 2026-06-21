# F07 Training Enrichment — Local Style Overrides
# Sources shared/style.R then adds F07-specific palettes and helpers.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

# ── Significance classification (Training concordance) ──────────────────────

SIG_COLORS_F7 <- c(

  "Interaction"    = "#9B7FBF",
  "Sig Both"       = "#2E7D32",
  "Sig HR only"    = "#2166AC",
  "Sig LR only"    = "#B2182B",
  "NS"             = "grey70"
)

SIG_LABEL_FILL_F7 <- c(
  "Interaction"    = scales::alpha("#9B7FBF", 0.75),
  "Sig Both"       = scales::alpha("#2E7D32", 0.75),
  "Sig HR only"    = scales::alpha("#2166AC", 0.75),
  "Sig LR only"    = scales::alpha("#B2182B", 0.75),
  "NS"             = scales::alpha("grey70",  0.75)
)

SIG_LABEL_TEXT_F7 <- c(
  "Interaction"    = "white",
  "Sig Both"       = "white",
  "Sig HR only"    = "white",
  "Sig LR only"    = "white",
  "NS"             = "white"
)

classify_proteins_f7 <- function(pi_HR, pi_LR, pi_int, threshold = 0.05) {
  dplyr::case_when(
    pi_int < threshold                      ~ "Interaction",
    pi_HR < threshold & pi_LR < threshold   ~ "Sig Both",
    pi_HR < threshold                       ~ "Sig HR only",
    pi_LR < threshold                       ~ "Sig LR only",
    TRUE                                    ~ "NS"
  ) |>
    factor(levels = c("Interaction", "Sig Both",
                       "Sig HR only", "Sig LR only", "NS"))
}

# ── Display labels for ORA bar pathway names ─────────────────────────────────

DISPLAY_LABELS_F07 <- c(
  "Formation Of The Dystrophin Glycoprotein Complex Dgc" = "Dystrophin Complex",
  "Striated Muscle Contraction"        = "Striated Muscle Contraction",
  "Reference Translation Initiation"   = "Translation Initiation",
  "Respiratory Electron Transport"     = "Respiratory ETC",
  "Complex I Biogenesis"               = "Complex I Biogenesis",
  "Epithelial Mesenchymal Transition"   = "EMT",
  "ATP Synthesis Coupled Electron Transport" = "ATP Synthesis e- Transport"
)

# ── Sigmoid ribbon (for sankey in panel B) ──────────────────────────────────

make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble::tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id
  )
}
