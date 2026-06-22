# panel_A_volcanoes.R — HRvLR volcano-ring panels (main + supplementary)
# Main figure: 7 state/within-group contrasts
# Supplement:  2 interaction contrasts
# Adapted from the YvO F04 per-panel export contract.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggforce)
})

here::i_am(".here")
source(here::here("04_Figures/shared/pathway_utils.R"))
source(here::here("04_Figures/F04/a_script/style.R"))
source(here::here("04_Figures/F04/a_script/volcano_ring.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────
RPT_MAIN_PDF <- here::here("04_Figures/F04/b_reports/main/pdf")
RPT_MAIN_PNG <- here::here("04_Figures/F04/b_reports/main/png")
RPT_SUPP_PDF <- here::here("04_Figures/F04/b_reports/supp/pdf")
RPT_SUPP_PNG <- here::here("04_Figures/F04/b_reports/supp/png")
DAT <- here::here("04_Figures/F04/c_data")
for (d in c(RPT_MAIN_PDF, RPT_MAIN_PNG, RPT_SUPP_PDF, RPT_SUPP_PNG, DAT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Load data once ────────────────────────────────────────────────────────────
dep_df <- read_csv(
  here::here("03_DEP/a_non_imputed/c_data/03_combined_results.csv"),
  show_col_types = FALSE
)

# ── Panel dimensions ─────────────────────────────────────────────────────────
PW <- 190; PH <- 185

# ── Contrast config ──────────────────────────────────────────────────────────
main_specs <- list(
  list(panel = "panel_A", contrast = "Baseline_HRvLR",
       title = "Baseline State: HR vs LR",
       subtitle = "HR - LR at Baseline (T1)"),
  list(panel = "panel_B", contrast = "Trained_HRvLR",
       title = "Trained State: HR vs LR",
       subtitle = "HR - LR after training (T2)"),
  list(panel = "panel_C", contrast = "Acute_HRvLR",
       title = "Acute State: HR vs LR",
       subtitle = "HR - LR after acute bout (T3)"),
  list(panel = "panel_D", contrast = "Training_HR",
       title = "Training Response - HR",
       subtitle = "HR: T2 - T1 (72h post-training)"),
  list(panel = "panel_E", contrast = "Training_LR",
       title = "Training Response - LR",
       subtitle = "LR: T2 - T1 (72h post-training)"),
  list(panel = "panel_F", contrast = "Acute_HR",
       title = "Acute Response - HR",
       subtitle = "HR: T3 - T2 (1h acute bout)"),
  list(panel = "panel_G", contrast = "Acute_LR",
       title = "Acute Response - LR",
       subtitle = "LR: T3 - T2 (1h acute bout)")
)

supp_specs <- list(
  list(panel = "panel_A", contrast = "Training_Interaction",
       title = "Training Interaction (HR - LR)",
       subtitle = "(HR_T2 - HR_T1) - (LR_T2 - LR_T1)"),
  list(panel = "panel_B", contrast = "Acute_Interaction",
       title = "Acute Interaction (HR - LR)",
       subtitle = "(HR_T3 - HR_T2) - (LR_T3 - LR_T2)")
)

all_specs <- c(main_specs, supp_specs)

# ── Build local fGSEA cache once from current DEP results ────────────────────
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500)

fgsea_df <- purrr::map_dfr(all_specs, function(spec) {
  contrast <- spec$contrast
  stat_col <- paste0("t_", contrast)
  ranks <- setNames(dep_df[[stat_col]], dep_df$gene)
  ranks <- ranks[!is.na(ranks) & !duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)

  run_fgsea_deduplicated(
    ranks = ranks,
    pathways = pw_collection,
    jaccard_cutoff = 0.5,
    nperm = 10000,
    min_size = 15,
    max_size = 500
  ) |>
    mutate(
      contrast = contrast,
      leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")
    )
})

write_csv(fgsea_df, file.path(DAT, "01_fgsea_cache.csv"))

build_and_save_panel <- function(spec, report_pdf_dir, report_png_dir,
                                 data_root, file_prefix) {
  contrast <- spec$contrast
  panel_dat <- file.path(data_root, spec$panel)
  dir.create(panel_dat, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("  Building %s panel %s (%s)", file_prefix, spec$panel, contrast))

  top_terms <- select_ring_terms(fgsea_df, contrast, n_each = 6)
  ring_data <- build_ring_with_gaps(top_terms, contrast, fgsea_df, n_each = 6)

  p <- make_volcano_ring(
    de_df = dep_df,
    go_df = fgsea_df,
    contrast = contrast,
    ring_data_override = ring_data,
    contrast_title = spec$title,
    contrast_subtitle = spec$subtitle,
    title_size = scale_text(BASE_TAG, PW),
    label_size = scale_text(BASE_PATHWAY, PW),
    count_label_size = scale_text(BASE_COUNT, PW),
    point_size = 1.2,
    point_alpha = 0.55
  )

  save_panel(
    p,
    file.path(report_png_dir, paste0(file_prefix, "_", spec$panel, "_", contrast)),
    PW, PH
  )
  file.copy(
    file.path(report_png_dir, paste0(file_prefix, "_", spec$panel, "_", contrast, ".pdf")),
    file.path(report_pdf_dir, paste0(file_prefix, "_", spec$panel, "_", contrast, ".pdf")),
    overwrite = TRUE
  )

  ring_out <- attr(p, "ring_data")
  if (!is.null(ring_out) && nrow(ring_out) > 0) {
    write_csv(
      ring_out |> dplyr::select(-gene_list),
      file.path(panel_dat, "ring_terms.csv")
    )
  }

  dep_df |>
    transmute(
      gene,
      logFC = round(.data[[paste0("logFC_", contrast)]], 4),
      neg_log10_p = round(-log10(.data[[paste0("P.Value_", contrast)]]), 4),
      pi_score = round(.data[[paste0("pi_score_", contrast)]], 6),
      adj_pval = round(.data[[paste0("adj.P.Val_", contrast)]], 6),
      direction = case_when(
        .data[[paste0("pi_score_", contrast)]] < 0.05 &
          .data[[paste0("logFC_", contrast)]] > 0 ~ "Up",
        .data[[paste0("pi_score_", contrast)]] < 0.05 &
          .data[[paste0("logFC_", contrast)]] < 0 ~ "Down",
        TRUE ~ "NS"
      )
    ) |>
    filter(!is.na(logFC), !is.na(neg_log10_p)) |>
    arrange(pi_score) |>
    write_csv(file.path(panel_dat, paste0("volcano_", contrast, ".csv")))
}

walk(
  main_specs,
  build_and_save_panel,
  report_pdf_dir = RPT_MAIN_PDF,
  report_png_dir = RPT_MAIN_PNG,
  data_root = file.path(DAT, "main"),
  file_prefix = "MAIN"
)

walk(
  supp_specs,
  build_and_save_panel,
  report_pdf_dir = RPT_SUPP_PDF,
  report_png_dir = RPT_SUPP_PNG,
  data_root = file.path(DAT, "supp"),
  file_prefix = "SUPP"
)

message("All F04 volcano-ring panels complete.")
