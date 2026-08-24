# The candidate labels the stage-04 sweep fits. Six here; stage 03 appends a
# seventh if the proteome turns out to carry its own split.
#
# Each label is a median cut into hi and lo. Dichotomising a continuous outcome
# costs power (Cohen 1983; Royston, Altman & Sauerbrei 2006) and the continuous
# form of this question was already answered in V1 across 54 association cells
# with nothing promoted. These labels exist to put every candidate on the same
# footing as the given HR/LR label, which is itself a median cut and cannot be
# compared against a continuous fit.
#
# The internal flag records whether an outcome is an ingredient of
# comp_hypertrophy. An internal label cannot corroborate the composite that
# contains it, so its sweep result is descriptive only. Three labels are
# genuinely external: both 1RM measures and volume load, the last of which
# shares r2 = 0.002 with the composite and describes what a subject did rather
# than what happened to them.

pacman::p_load(here, dplyr, tidyr, purrr, tibble, readr, openxlsx)

OUT_DIR <- here("03_Features", "01_Responsiveness", "c_data")

INTERNAL_R2 <- 0.20

pheno <- read_csv(here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
)
structure_tbl <- openxlsx::read.xlsx(
  file.path(OUT_DIR, "01_label_audit.xlsx"), "composite_structure"
)

median_split <- function(x, subject) {
  keep <- !is.na(x)
  tibble(
    subject = subject[keep],
    level = ifelse(x[keep] > stats::median(x[keep]), "hi", "lo")
  )
}

trait_labels <- map_dfr(structure_tbl$trait, function(v) {
  median_split(pheno[[v]], pheno$subject) |>
    mutate(label = sub("^d_", "", v), .before = 1)
})

given <- tibble(
  label = "given",
  subject = pheno$subject,
  level = ifelse(pheno$group_arm == "HR", "hi", "lo")
)

candidate_labels <- bind_rows(given, trait_labels)

# Several labels land on the same split: the three fibre-area measures
# correlate above 0.9 and their median cuts are identical. Fitting each would
# report one test three times and inflate the sweep's apparent size, so each
# partition is named once and the labels sharing it are recorded against it. A
# split and its complement are the same partition, hence the canonical form.
partition_of <- function(level, subject) {
  ord <- order(subject)
  s <- paste(level[ord], collapse = "")
  min(s, chartr("hilo", "lohi", s))
}

partitions <- candidate_labels |>
  summarise(key = partition_of(.data$level, .data$subject), .by = label) |>
  mutate(partition = as.integer(factor(.data$key, levels = unique(.data$key))))

labels_meta <- candidate_labels |>
  summarise(
    n = dplyr::n(),
    n_hi = sum(level == "hi"),
    n_lo = sum(level == "lo"),
    .by = label
  ) |>
  left_join(
    structure_tbl |>
      transmute(
        label = sub("^d_", "", trait),
        r_with_composite = r, r2_alone
      ),
    by = "label"
  ) |>
  mutate(
    r_with_composite = dplyr::coalesce(r_with_composite, 1),
    r2_alone = dplyr::coalesce(r2_alone, 1),
    internal = r2_alone >= INTERNAL_R2,
    agreement_with_given = map_dbl(label, function(l) {
      j <- inner_join(
        filter(candidate_labels, label == l),
        filter(candidate_labels, label == "given"),
        by = "subject", suffix = c("", "_given")
      )
      mclust::adjustedRandIndex(j$level, j$level_given)
    })
  ) |>
  left_join(
    partitions |> dplyr::select(label, partition),
    by = "label"
  ) |>
  mutate(
    shares_partition_with = purrr::map_chr(seq_along(.data$label), function(i) {
      other <- .data$label[.data$partition == .data$partition[i]]
      paste(setdiff(other, .data$label[i]), collapse = ", ")
    })
  ) |>
  arrange(.data$partition, desc(.data$r2_alone))

write_csv(candidate_labels, file.path(OUT_DIR, "02_candidate_labels.csv"))
write.xlsx(
  list(candidate_labels = candidate_labels, labels_meta = labels_meta),
  file.path(OUT_DIR, "02_candidate_labels.xlsx")
)

print(as.data.frame(labels_meta), digits = 3)
message(
  "\n", nrow(labels_meta), " labels resolve to ",
  dplyr::n_distinct(labels_meta$partition), " distinct partitions"
)
