# Pathway level, two engines answering two different questions.
#
# fgsea reads the whole ranked list, so it can find a coordinated shift across a
# set whose members individually clear nothing. That is why this stage runs
# despite stage 04's three BH survivors failing their permutation check: the
# survivors are not what fgsea is built on.
#
# singscore scores each sample independently and is then fitted through the same
# estimator the proteins used, so a pathway logFC and a protein logFC mean the
# same thing and sit in the same table. Its scores are computed once upstream
# and cached; recomputing them per label would be the same numbers each time.
#
# No fry. BH within each label, contrast and engine, never pooled across the
# twelve cells.

pacman::p_load(
  here, dplyr, tidyr, purrr, tibble, readr, fgsea, openxlsx
)

source(here("functions", "label_contrasts.R"))
source(here("functions", "shared_pathway_utils.R"))

OUT_DIR <- here("03_Features", "05_Pathways", "c_data")
BH_ALPHA <- 0.05

sweep <- read_csv(
  here("03_Features", "04_Proteins", "c_data", "01_sweep_full.csv"),
  show_col_types = FALSE
)
labels_long <- read_csv(
  here("03_Features", "01_Responsiveness", "c_data", "02_candidate_labels.csv"),
  show_col_types = FALSE
)
annotation <- read_csv(
  here("02_Normalization", "c_data", "normalized.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(feature = uniprot_id, gene)

pathways <- build_pathway_collection()

# fgsea ranks by gene, and several uniprot ids can carry one symbol. The most
# extreme statistic wins rather than an average, which would pull a real shift
# toward zero whenever one isoform is flat.
rank_vector <- function(df) {
  df |>
    left_join(annotation, by = "feature") |>
    filter(!is.na(.data$gene), .data$gene != "", !is.na(.data$t)) |>
    slice_max(abs(.data$t), n = 1, by = "gene", with_ties = FALSE) |>
    (\(d) setNames(d$t, d$gene))()
}

cells <- sweep |> distinct(label, contrast)

set.seed(42)
fgsea_res <- pmap_dfr(cells, function(label, contrast) {
  ranks <- rank_vector(filter(
    sweep, .data$label == !!label, .data$contrast == !!contrast
  ))
  run_fgsea(ranks, pathways) |>
    mutate(label = !!label, contrast = !!contrast, .before = 1)
})

singscore_res <- map_dfr(unique(labels_long$label), function(l) {
  lab <- filter(labels_long, .data$label == l)
  vec <- setNames(lab$level, lab$subject)
  fit_label_contrasts(pathway_matrix(), vec) |>
    mutate(label = l, .before = 1)
})

fgsea_summary <- fgsea_res |>
  summarise(
    n_sets = dplyr::n(),
    n_bh = sum(.data$padj < BH_ALPHA, na.rm = TRUE),
    min_padj = min(.data$padj, na.rm = TRUE),
    .by = c(label, contrast)
  )

singscore_summary <- singscore_res |>
  summarise(
    n_sets = dplyr::n(),
    n_bh = sum(.data$bh < BH_ALPHA, na.rm = TRUE),
    min_bh = min(.data$bh, na.rm = TRUE),
    .by = c(label, contrast)
  )

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(fgsea_res, file.path(OUT_DIR, "01_fgsea_full.csv"))
write.xlsx(
  list(
    fgsea_summary = fgsea_summary,
    fgsea_hits = filter(fgsea_res, .data$padj < BH_ALPHA) |> arrange(padj),
    singscore_summary = singscore_summary,
    singscore_hits = filter(singscore_res, .data$bh < BH_ALPHA) |> arrange(bh)
  ),
  file.path(OUT_DIR, "01_pathways.xlsx")
)

print(as.data.frame(fgsea_summary), digits = 3)
message(
  "\nfgsea: ", sum(fgsea_summary$n_bh), " set-contrast hits at BH < ", BH_ALPHA,
  " | singscore: ", sum(singscore_summary$n_bh)
)
