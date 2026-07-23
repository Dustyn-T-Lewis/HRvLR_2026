#!/usr/bin/env Rscript
# Builds the red-cell proteome reference consumed by Stage 01's contamination filter. The five
# source datasets ship in RBC_proteins.xlsx; this writes the collapsed lookup and nothing else,
# so it belongs at stage 00, upstream of filtering.

pacman::p_load(here, readxl, readr, dplyr, tidyr, stringr, purrr)
source(here("00_input", "a_script", "build_blood_list.R"))

rbc_reference <- build_rbc_reference(here("00_input", "RBC_proteins.xlsx"))

write_tsv(rbc_reference, here("00_input", "RBC_proteome_reference.tsv"))

cat(sprintf(
  "RBC_proteome_reference.tsv: %d proteins (%d with accession, %d with gene); %d backed by >=2 sources\n",
  nrow(rbc_reference), sum(!is.na(rbc_reference$acc)), sum(!is.na(rbc_reference$gene)),
  sum(rbc_reference$n_sources >= 2)
))
