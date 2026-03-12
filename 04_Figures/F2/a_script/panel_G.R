# panel_G.R — Acute_Interaction Volcano + Ring
# Contrast: Acute_Interaction

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggforce)
})

here::i_am(".here")
source(here::here("04_Figures/F2/a_script/style_f2.R"))
source(here::here("04_Figures/F2/a_script/volcano_ring.R"))

CONTRAST <- "Acute_Interaction"
PW <- 190; PH <- 185
RPT <- here::here("04_Figures/F2/b_reports")
DAT <- here::here("04_Figures/F2/c_data")
dir.create(file.path(DAT, "panel_G"), recursive = TRUE, showWarnings = FALSE)

dep_df   <- read_csv(here::here("03_DEP/c_data/03_combined_results.csv"),
                     show_col_types = FALSE)
fgsea_df <- read_csv(here::here("04_Figures/F1/c_data/06_panel_F_fgsea_results.csv"),
                     show_col_types = FALSE)

top_terms <- select_ring_terms(fgsea_df, CONTRAST, n_each = 6)
ring_data  <- build_ring_with_gaps(top_terms, CONTRAST, fgsea_df, n_each = 6)

p <- make_volcano_ring(
  de_df              = dep_df,
  go_df              = fgsea_df,
  contrast           = CONTRAST,
  ring_data_override = ring_data,
  contrast_title     = "Acute Interaction (HR \u2212 LR)",
  contrast_subtitle  = "Acute: HR \u2212 LR response",
  title_size         = scale_text(BASE_TAG, PW),
  label_size         = scale_text(BASE_PATHWAY, PW),
  count_label_size   = scale_text(BASE_COUNT, PW),
  point_size         = 1.2,
  point_alpha        = 0.55
)

save_panel(p, file.path(RPT, "panel_G_Acute_Interaction"), PW, PH)

# Save ring terms used
ring_out <- attr(p, "ring_data")
if (!is.null(ring_out) && nrow(ring_out) > 0)
  write_csv(ring_out |> dplyr::select(-gene_list),
            file.path(DAT, "panel_G", "ring_terms.csv"))

# Save volcano data
dep_df |>
  dplyr::transmute(
    gene,
    logFC         = round(.data[[paste0("logFC_",        CONTRAST)]], 4),
    neg_log10_p   = round(-log10(.data[[paste0("P.Value_", CONTRAST)]]), 4),
    pi_score      = round(.data[[paste0("pi_score_",     CONTRAST)]], 6),
    adj_pval      = round(.data[[paste0("adj.P.Val_",    CONTRAST)]], 6),
    direction     = dplyr::case_when(
      .data[[paste0("pi_score_", CONTRAST)]] < 0.05 &
        .data[[paste0("logFC_",    CONTRAST)]] > 0 ~ "Up",
      .data[[paste0("pi_score_", CONTRAST)]] < 0.05 &
        .data[[paste0("logFC_",    CONTRAST)]] < 0 ~ "Down",
      TRUE ~ "NS")
  ) |>
  dplyr::filter(!is.na(logFC), !is.na(neg_log10_p)) |>
  dplyr::arrange(pi_score) |>
  write_csv(file.path(DAT, "panel_G", paste0("volcano_Acute_Interaction.csv")))

message("Panel G done: Acute_Interaction")
