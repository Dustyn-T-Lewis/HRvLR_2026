# Pathways from each trait's protein ranking, every subject pooled. Proteins are ranked by the limma
# t that 01_Associate/01_Proteins saved for each pairing and trait; fgsea and cameraPR run on that
# ranking. fgsea's p and BH are its own; cameraPR is reported beside it, with camera_agrees marking
# sets fgsea calls at FDR < 0.05 that cameraPR reaches at p < 0.05 in the same direction.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})
source(here("R", "helpers.R"))

stage <- here("04_Continuous_no_group", "02_Pathways")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  associate = "04_Continuous_no_group/01_Associate/01_Proteins/c_data/associate.rds",
  gene_sets = "01_Trunk/04_Gene_Sets/c_data/gene_sets.rds"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

gs <- readRDS(paths[["gene_sets"]])
catalog <- filter(gs$set_catalog, qualifies)
selected <- filter(gs$protein_map, selected)
gene_of <- set_names(selected$gene, selected$uniprot_id)
protein_t <- readRDS(paths[["associate"]])$change_score |>
  filter(feature %in% selected$uniprot_id) |>
  select(window, outcome, feature, t)

rank_tests <- function(d) {
  stats <- set_names(d$t, gene_of[d$feature])
  set.seed(1)
  fg <- fgsea::fgsea(gs$sets, stats, minSize = 15, maxSize = 500) |>
    as_tibble() |>
    transmute(set_id = pathway, size, nes = NES, p = pval, fdr = padj, leading_edge = leadingEdge)
  index <- map(gs$sets, \(genes) which(names(stats) %in% genes))
  camera <- cameraPR(stats, index[lengths(index) >= 15]) |>
    as_tibble(rownames = "set_id") |>
    transmute(set_id, camera_direction = Direction, p_camera = PValue)
  left_join(fg, camera, by = "set_id") |>
    mutate(camera_agrees = fdr < 0.05 & p_camera < 0.05 & sign(nes) == if_else(
      camera_direction == "Up", 1, -1
    ))
}
pathways <- protein_t |>
  nest(.by = c(window, outcome)) |>
  mutate(tests = map(data, rank_tests), .keep = "unused") |>
  unnest(tests) |>
  left_join(select(catalog, set_id, collection, pathway), by = "set_id") |>
  relocate(collection, pathway, .after = set_id) |>
  arrange(window, outcome, p)

counts <- pathways |>
  summarise(
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05, na.rm = TRUE),
    n_up = sum(fdr < 0.05 & nes > 0, na.rm = TRUE),
    n_down = sum(fdr < 0.05 & nes < 0, na.rm = TRUE),
    n_camera_agrees = sum(camera_agrees, na.rm = TRUE),
    .by = c(window, outcome)
  )

top_sets <- pathways |>
  filter(fdr < 0.05) |>
  slice_min(p, n = 10, by = c(window, outcome), with_ties = FALSE) |>
  mutate(label = stringr::str_trunc(pathway, 60))
pdf(file.path(figures, "02_pathways_figures.pdf"), width = 11, height = 8.5)
print(
  ggplot(
    tidyr::pivot_longer(counts, c(n_up, n_down), names_to = "direction"),
    aes(value, outcome, fill = direction)
  ) +
    geom_col() +
    facet_wrap(~window, nrow = 1) +
    scale_fill_manual(
      values = c(n_up = "#4D9221", n_down = "#C51B7D"), labels = c(n_up = "up", n_down = "down")
    ) +
    labs(
      title = "Sets at fgsea FDR < 0.05 per trait ranking",
      subtitle = "Proteins ranked by limma t against each trait; BH within window and trait",
      x = "Sets", y = NULL, fill = NULL
    ) +
    theme_bw(base_size = 9)
)
walk(split(top_sets, top_sets$window), \(d) {
  print(
    ggplot(d, aes(nes, reorder(paste(label, outcome, sep = "___"), nes))) +
      geom_col(aes(fill = nes > 0), show.legend = FALSE) +
      scale_y_discrete(labels = \(x) sub("___.*$", "", x)) +
      facet_wrap(~outcome, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c(`TRUE` = "#4D9221", `FALSE` = "#C51B7D")) +
      labs(
        title = sprintf("Top sets per trait, %s window", d$window[1]),
        subtitle = "Up to ten sets at fgsea FDR < 0.05 per trait, smallest p first",
        x = "NES", y = NULL
      ) +
      theme_bw(base_size = 7)
  )
})
invisible(dev.off())

write_workbook(
  list(counts = counts, pathways = select(pathways, -leading_edge)),
  c(
    "Per window and trait: sets tested, at p < 0.05, at FDR < 0.05 (up, down), cameraPR agreeing.",
    "Per window, trait and set: fgsea NES, p, BH; cameraPR direction and p; camera_agrees."
  ),
  file.path(out, "02_pathways.xlsx"), inputs, paths
)
saveRDS(pathways, file.path(out, "pathways.rds"), compress = "xz")
