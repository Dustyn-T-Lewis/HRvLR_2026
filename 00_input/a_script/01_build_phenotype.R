#!/usr/bin/env Rscript
# Writes the canonical per-subject trait table. Its only input is the sample sheet, so it
# belongs at stage 00, upstream of every consumer (F01, the F04 association figures, F05). Before this
# existed, 00_input/c_data/phenotype.csv was an orphan: `modules` read it and nothing built it.

pacman::p_load(here, readr)
source(here("00_input", "a_script", "build_phenotype.R"))

phenotype <- build_phenotype_table(here("00_input", "HRvLR_meta.csv"))

out <- here("00_input", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv(phenotype, file.path(out, "phenotype.csv"))

cat(sprintf(
  "phenotype.csv: %d subjects x %d traits\n",
  nrow(phenotype), ncol(phenotype) - 2
))
