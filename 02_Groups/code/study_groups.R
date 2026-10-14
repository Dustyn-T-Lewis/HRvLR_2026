#' ---
#' title: "Study groups: HR and LR"
#' params:
#'   grouping: HRvLR
#' format:
#'   typst:
#'     toc: true
#'     margin:
#'       x: 1.6cm
#' execute:
#'   echo: true
#'   warning: false
#' ---
#'
#' The study's own grouping. 25 people trained; each got a composite hypertrophy score; the 8 with
#' the highest scores were profiled as higher responders (HR) and the 8 with the lowest as lower
#' responders (LR). The middle 9 were not profiled. This is an extreme-groups design: the two groups
#' sit far apart by construction, which makes a difference easier to see and makes any correlation
#' computed across these 16 people stronger than it would be across all 25.
#'
#' This step checks the labels against the phenotype and writes them in the same form as the
#' clustering groupings, so every comparison downstream reads one kind of groups file.
#'
#' ## Packages

#| label: packages
#| message: false
library(here)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))

#' ## Inputs
#'
#' `phenotype.csv` holds one row per person: muscle and fibre cross-sectional area (mCSA, fCSA) and
#' 1RM strength before and after training, the composite score and the training volume.
#' `metadata.csv` maps each biopsy run to its person.

#| label: inputs
out <- here("02_Groups", params$grouping, "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  phenotype = "00_Setup/Input/phenotype.csv",
  samples = "00_Setup/Input/metadata.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
samples <- readr::read_csv(paths[["samples"]], show_col_types = FALSE)
count(phenotype, group)

#' ## The composite score
#'
#' The manuscript defines the composite as (2 × %ΔmCSA + %ΔfCSA) / 3, weighting whole-muscle growth
#' twice as heavily as mixed fibre growth. Recomputed from the pre and post columns, it must match
#' the `comp_hypertrophy` column the source supplied, up to the source's rounding to two decimals.

#| label: composite
composite <- phenotype |>
  transmute(
    subject, group,
    pct_mcsa = 100 * (mcsa_post_cm2 / mcsa_pre_cm2 - 1),
    pct_fcsa = 100 * (fcsa_mixed_post_um2 / fcsa_mixed_pre_um2 - 1),
    recomputed = (2 * pct_mcsa + pct_fcsa) / 3,
    comp_hypertrophy
  )
stopifnot(all(abs(composite$recomputed - composite$comp_hypertrophy) < 0.01 + 1e-9))
composite |>
  summarise(
    people = n(), mean = mean(comp_hypertrophy), sd = sd(comp_hypertrophy),
    lowest = min(comp_hypertrophy), highest = max(comp_hypertrophy), .by = group
  )


#| label: fig-composite
#| fig-cap: "Composite score per person. The 9 unprofiled people sat in the gap."
#| fig-width: 5
#| fig-height: 3.5
ggplot(composite, aes(group, comp_hypertrophy, colour = group)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_point(position = position_jitter(width = 0.1, seed = 1), size = 2) +
  scale_colour_manual(values = c(HR = "#2166AC", LR = "#B2182B"), guide = "none") +
  labs(x = NULL, y = "composite hypertrophy (%)") +
  theme_minimal(base_size = 10)

#' ## Labels against the sample sheet
#'
#' Every person in the phenotype table must carry the same label in the sample sheet.

#| label: labels
labels_in_samples <- distinct(samples, subject, group)
group_in_phenotype <- with(phenotype, set_names(group, subject))
stopifnot(
  !anyDuplicated(labels_in_samples$subject),
  setequal(labels_in_samples$subject, phenotype$subject),
  all(labels_in_samples$group == group_in_phenotype[labels_in_samples$subject])
)
group_of <- with(labels_in_samples, set_names(group, subject))
table(group_of)

#' ## Summary
#'
#' `r sum(group_of == "HR")` HR and `r sum(group_of == "LR")` LR people. HR's composite runs from
#' `r round(min(composite$comp_hypertrophy[composite$group == "HR"]), 2)` to
#' `r round(max(composite$comp_hypertrophy[composite$group == "HR"]), 2)`%, LR's from
#' `r round(min(composite$comp_hypertrophy[composite$group == "LR"]), 2)` to
#' `r round(max(composite$comp_hypertrophy[composite$group == "LR"]), 2)`%. Because the composite
#' defines the groups, a difference in muscle or fibre growth between HR and LR is built in and is
#' not a finding; `comp_hypertrophy` stays out of the group-by-trait tests.
#'
#' ## Save
#'
#' `groups.rds` holds who is in each group, the group order (the first group comes first in every
#' contrast) and the traits that defined the groups.

#| label: save
saveRDS(
  list(group_of = group_of, groups = c("HR", "LR"), defining_traits = "comp_hypertrophy"),
  file.path(out, "groups.rds")
)
write_workbook(
  list(composite = composite),
  "Per person: % change in mCSA and mixed fCSA, the recomputed composite and the source's.",
  file.path(out, "study_groups.xlsx"), inputs, paths
)
