# Writes the README of every analysis step from its saved outputs: what it read (its workbooks'
# input_manifest), what it wrote (the files on disk), how to run it, its key result (computed from
# its workbook) and its sheets (the workbook's read_me). The question, method and cautions come
# from docs/readme_text.R. Folder and figure READMEs are written by hand. Run after run_all.R.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(purrr)
  library(tibble)
})
source(here("docs", "readme_text.R"))

run_all <- parse(here("run_all.R"))
stages_at <- which(map_lgl(as.list(run_all), \(e) identical(e[[2]], quote(stages))))
stages <- eval(run_all[[stages_at]][[3]])
inputs_used <- c("T1", "T2", "T3", "training", "acute")

rd <- function(file, sheet) readxl::read_excel(file, sheet)
workbook <- function(stage) list.files(here(stage, "c_data"), "[.]xlsx$", full.names = TRUE)[1]
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
p_text <- function(p) if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)
range_by <- function(d, n, label) {
  lo <- d[which.min(d[[n]]), ]
  hi <- d[which.max(d[[n]]), ]
  sprintf("%s (%s) to %s (%s)", fmt(lo[[n]]), lo[[label]], fmt(hi[[n]]), hi[[label]])
}

# Key results, one function per kind of step. Each returns markdown bullet lines.
key_filtering <- function(stage) {
  f <- workbook(stage)
  log <- rd(f, "filter_log")
  runs <- rd(f, "outlier_diagnostics")
  dropped <- runs$sample_id[runs$consensus_outlier]
  c(
    sprintf(
      paste(
        "- %s proteins in the raw report; %d removed as contaminants, %d by the detection filter;",
        "%s kept."
      ),
      fmt(log$n_after[1]), sum(log$n_removed[grepl("^remove", log$step)]),
      log$n_removed[log$step == "detection"], fmt(tail(log$n_after, 1))
    ),
    sprintf(
      "- %d of %d runs fail three or more of the four checks and are dropped: %s.",
      length(dropped), nrow(runs), paste(dropped, collapse = ", ")
    )
  )
}
key_normalisation <- function(stage) {
  f <- workbook(stage)
  eta <- rd(f, "eta_squared")
  pcs <- rd(f, "components")
  c(
    sprintf(
      paste(
        "- %s proteins normalised; the six group-by-time cells explain a median %.0f%% of a",
        "protein's variance (eta squared)."
      ),
      fmt(nrow(eta)), 100 * median(eta$eta2, na.rm = TRUE)
    ),
    sprintf(
      "- PC1 and PC2 carry %.1f%% and %.1f%% of the variance.",
      100 * pcs$variance_explained[1], 100 * pcs$variance_explained[2]
    )
  )
}
key_quantification <- function(stage) {
  f <- workbook(stage)
  curve <- rd(f, "curve")
  proteins <- rd(f, "proteins")
  c(
    sprintf(
      paste(
        "- %.1f%% of values are missing; the detection curve has intercept %.2f and slope %.2f",
        "(logit scale)."
      ),
      curve$pct_missing[1], curve$estimate[curve$term == "beta0"],
      curve$estimate[curve$term == "beta1"]
    ),
    sprintf(
      paste(
        "- %s proteins get a value in all %d biopsies; %s were measured in every one, and the",
        "median protein in %d."
      ),
      fmt(nrow(proteins)), max(proteins$n_observed), fmt(sum(proteins$prop_observed == 1)),
      as.integer(median(proteins$n_observed))
    )
  )
}
key_gene_sets <- function(stage) {
  collections <- rd(workbook(stage), "collection_summary")
  mapping <- rd(workbook(stage), "mapping_summary")
  c(
    sprintf(
      "- %s sets qualify for testing: %s.", fmt(sum(collections$n_qualifying)),
      paste(
        sprintf(
          "%s %s of %s", collections$collection, fmt(collections$n_qualifying),
          fmt(collections$n_in_source)
        ),
        collapse = ", "
      )
    ),
    if (identical(mapping$mapping_status, "representative")) {
      sprintf("- All %s proteins map to one distinct gene symbol each.", fmt(mapping$n_proteins))
    } else {
      sprintf(
        "- Mapping: %s.",
        paste(sprintf("%s `%s`", fmt(mapping$n_proteins), mapping$mapping_status), collapse = ", ")
      )
    }
  )
}
key_set_scores <- function(stage) {
  f <- workbook(stage)
  scores <- rd(f, "score_summary")
  structure <- rd(f, "structure_check")
  c(
    sprintf("- %s sets scored in %d samples.", fmt(scores$n_sets), scores$n_samples),
    sprintf(
      "- Subject identity explains %.0f%% of PC1 and %.0f%% of PC2 of the scores.",
      100 * structure$subject_share[1], 100 * structure$subject_share[2]
    )
  )
}
key_module_building <- function(stage) {
  f <- workbook(stage)
  modules <- filter(rd(f, "module_summary"), module != "grey")
  membership <- rd(f, "membership")
  power <- filter(rd(f, "soft_threshold"), chosen)
  c(
    sprintf(
      "- %d modules of %d to %d proteins; %d proteins fit no module (`grey`).",
      nrow(modules), min(modules$n_proteins), max(modules$n_proteins),
      sum(membership$module == "grey")
    ),
    sprintf("- Soft power %d (scale-free fit R² %.2f).", power$Power, power$SFT.R.sq)
  )
}
key_characterisation <- function(stage) {
  labels <- rd(workbook(stage), "labels")
  strongest <- slice_min(labels, top_set_fdr, n = 3)
  c(
    sprintf(
      paste(
        "- %d of %d modules have a gene set over-represented at FDR < 0.05; %d have more STRING",
        "edges than expected at FDR < 0.05."
      ),
      sum(labels$top_set_fdr < 0.05), nrow(labels), sum(labels$string_fdr < 0.05)
    ),
    sprintf(
      "- Strongest: %s.",
      paste(
        sprintf(
          "`%s`, %s (FDR %.1g)", strongest$module, strongest$top_set, strongest$top_set_fdr
        ),
        collapse = "; "
      )
    )
  )
}
key_design <- function(stage) {
  meta <- readRDS(here(stage, "c_data", "design.rds"))$dal$metadata
  per_arm <- meta |>
    summarise(people = n_distinct(subject), biopsies = n(), .by = arm) |>
    arrange(arm)
  correlation <- rd(workbook(stage), "correlation_strata")$correlation
  c(
    sprintf(
      "- %d biopsies from %d people: %s.", nrow(meta), n_distinct(meta$subject),
      paste(
        sprintf("%s %d people (%d biopsies)", per_arm$arm, per_arm$people, per_arm$biopsies),
        collapse = ", "
      )
    ),
    sprintf("- Within-subject correlation from limpa's fit: %.2f.", correlation)
  )
}
key_protein <- function(stage) {
  f <- workbook(stage)
  summary <- rd(f, "contrast_summary")
  permanova <- rd(f, "permanova")
  detectable <- rd(f, "detectable_effect")
  time <- filter(permanova, data == "proteome", term == "time", window == "all", arms == "both")
  group <- filter(permanova, data == "proteome", term == "arm", window == "subject mean")
  c(
    sprintf(
      paste(
        "- %s proteins tested in each of the nine contrasts; %s at FDR < 0.05 across all nine",
        "(lowest FDR %.3f)."
      ),
      fmt(summary$n_tested[1]), fmt(sum(summary$n_fdr_05)), min(summary$min_fdr)
    ),
    sprintf(
      "- Proteins at p < 0.05 run from %s; chance alone gives about %s per contrast.",
      range_by(summary, "n_nominal", "contrast"), fmt(round(summary$n_tested[1] * 0.05))
    ),
    sprintf(
      paste(
        "- PERMANOVA: time explains %.1f%% of the distance between biopsies (p %s); group, on one",
        "mean profile per person, %.1f%% (p %s)."
      ),
      100 * time$r2, p_text(time$p), 100 * group$r2, p_text(group$p)
    ),
    sprintf(
      paste(
        "- The Training_Interaction effect detectable at 80%% power is %.2f log2 units for the",
        "median protein."
      ),
      detectable$log2_fc[detectable$quantile == 0.5]
    )
  )
}
key_sets <- function(stage) {
  summary <- rd(workbook(stage), "set_summary")
  c(
    sprintf("- %s gene sets tested per contrast.", fmt(summary$n_tested_fgsea[1])),
    sprintf(
      "- fgsea calls %s sets at FDR < 0.05; the baseline contrast (the floor) calls %s.",
      range_by(summary, "n_fdr_fgsea", "contrast"),
      fmt(summary$n_fdr_fgsea[grepl("^Baseline", summary$contrast)])
    ),
    sprintf(
      "- fry calls %s sets at FDR < 0.05 across the nine contrasts.", fmt(sum(summary$n_fdr_fry))
    )
  )
}
key_modules <- function(stage) {
  f <- workbook(stage)
  eigengenes <- rd(f, "eigengene_tests")
  fry <- rd(f, "fry_tests")
  c(
    sprintf(
      "- %d modules by %d contrasts: %d eigengene tests at p < 0.05 and %d at FDR < 0.05.",
      n_distinct(eigengenes$module), n_distinct(eigengenes$contrast), sum(eigengenes$p < 0.05),
      sum(eigengenes$fdr < 0.05)
    ),
    sprintf(
      "- fry: %d module tests at p < 0.05 and %d at FDR < 0.05.", sum(fry$p < 0.05),
      sum(fry$fdr < 0.05)
    )
  )
}
key_preservation <- function(stage) {
  preservation <- rd(workbook(stage), "preservation")
  map_chr(split(preservation, preservation$reference), \(d) {
    sprintf(
      "- %s modules tested in %s: %d strongly preserved, %d weak to moderate, %d not preserved.",
      d$reference[1], d$test[1], sum(d$preservation == "strong"),
      sum(d$preservation == "weak to moderate"), sum(d$preservation == "none")
    )
  }) |>
    unname()
}
key_classification <- function(stage) {
  level <- sub(".*/0[1-3]_([A-Za-z]+).*", "\\1", sub("/(T[123]|training|acute)$", "", stage))
  noun <- c(Proteins = "proteins", Sets = "gene sets", Modules = "modules")[[level]]
  per_task <- rd(workbook(stage), "counts") |>
    summarise(
      n_nominal = sum(n_nominal), n_tested = sum(n_tested), n_fdr = sum(n_fdr_05), .by = task
    )
  c(
    sprintf(
      "- %s %s tested in each of %d tasks; %s at FDR < 0.05 in any task.",
      fmt(per_task$n_tested[1]), noun, nrow(per_task), fmt(sum(per_task$n_fdr))
    ),
    if (all(per_task$n_nominal == 0)) {
      sprintf("- No %s at p < 0.05 in any task.", sub("gene ", "", noun))
    } else {
      sprintf(
        "- %s at p < 0.05 run from %s.", tools::toTitleCase(noun),
        range_by(per_task, "n_nominal", "task")
      )
    }
  )
}
key_traits <- function(stage) {
  counts <- rd(workbook(stage), "counts")
  totals <- function(d) c(fmt(sum(d$n_tested)), fmt(sum(d$n_nominal)), fmt(sum(d$n_fdr_05)))
  difference <- totals(filter(counts, test == "group_difference"))
  within <- totals(filter(counts, test == "within_group"))
  c(
    sprintf(
      paste(
        "- Group difference (limma, feature ~ group * trait): %s tests, %s at p < 0.05, %s at",
        "FDR < 0.05."
      ),
      difference[1], difference[2], difference[3]
    ),
    sprintf(
      "- Within each group (Spearman): %s tests, %s at p < 0.05, %s at FDR < 0.05.",
      within[1], within[2], within[3]
    )
  )
}
key_clustering <- function(stage) {
  f <- workbook(stage)
  choice <- rd(f, "choice")
  stability <- rd(f, "stability")
  assigned <- filter(rd(f, "assignments"), !is.na(cluster))
  if (!"input" %in% names(choice)) {
    return(c(
      sprintf(
        "- %d people clustered; cluster sizes %s; bootstrap Jaccard %s.", choice$n_subjects,
        paste(count(assigned, cluster)$n, collapse = " and "),
        paste(sprintf("%.2f", filter(stability, k == choice$k_used)$jaccard), collapse = " and ")
      ),
      paste("-", choice$imposed_reason)
    ))
  }
  sizes <- count(assigned, input, cluster) |>
    summarise(sizes = paste(n, collapse = "/"), .by = input)
  jaccard <- filter(stability, k == 2) |>
    summarise(jaccard = paste(sprintf("%.2f", jaccard), collapse = "/"), .by = input)
  rows <- left_join(choice, sizes, by = "input") |> left_join(jaccard, by = "input")
  c(
    paste(
      "- The gap statistic selects K = 1 at every input; K = 2 is imposed so the two halves can",
      "be compared."
    ),
    "", "| Input | People | C1/C2 | Jaccard C1/C2 |", "|---|---|---|---|",
    sprintf("| %s | %d | %s | %s |", rows$input, rows$n_subjects, rows$sizes, rows$jaccard)
  )
}
key_cluster_phenotype <- function(stage) {
  tests <- rd(workbook(stage), "tests")
  sprintf(
    paste(
      "- %d cluster-against-trait tests (Wilcoxon) over the five inputs: %d at p < 0.05, %d at",
      "FDR < 0.05."
    ),
    nrow(tests), sum(tests$p < 0.05, na.rm = TRUE), sum(tests$fdr < 0.05, na.rm = TRUE)
  )
}
key_trait_associations <- function(stage) {
  names_of <- c(
    model = "Sample model", `change score` = "Change score", `two predictor` = "Two predictor",
    interaction = "Interaction"
  )
  rd(workbook(stage), "counts") |>
    mutate(kind = sub(":.*", "", analysis)) |>
    summarise(
      n_tested = sum(n_tested), n_nominal = sum(n_nominal), n_fdr = sum(n_fdr_05), .by = kind
    ) |>
    with(sprintf(
      "- %s: %s tests, %s at p < 0.05, %s at FDR < 0.05.", names_of[kind], fmt(n_tested),
      fmt(n_nominal), fmt(n_fdr)
    ))
}
key_trait_pathways <- function(stage) {
  counts <- rd(workbook(stage), "counts")
  by_window <- summarise(counts, n_fdr = sum(n_fdr_05), .by = window)
  top <- counts[which.max(counts$n_fdr_05), ]
  c(
    sprintf(
      "- %s sets ranked against each of %d trait and window pairs.", fmt(counts$n_tested[1]),
      nrow(counts)
    ),
    sprintf(
      "- Sets at FDR < 0.05, summed over traits: %s.",
      paste(sprintf("`%s` %s", by_window$window, fmt(by_window$n_fdr)), collapse = ", ")
    ),
    sprintf(
      "- The most for one trait: `%s` over `%s`, %s sets.", top$trait, top$window,
      fmt(top$n_fdr_05)
    )
  )
}
key_profiles <- function(stage) {
  f <- workbook(stage)
  groups <- rd(f, "profile_groups")
  trait_groups <- rd(f, "trait_groups")
  c(
    sprintf(
      "- %d profile groups hold %s proteins; the largest holds %d.", nrow(groups),
      fmt(sum(groups$n_proteins)), max(groups$n_proteins)
    ),
    sprintf(
      "- %d within-trait groups over %d trait and window pairs.", nrow(trait_groups),
      n_distinct(paste(trait_groups$window, trait_groups$trait))
    )
  )
}
key_of <- list(
  "01_Protein_Filtering$" = key_filtering, "02_Normalisation$" = key_normalisation,
  "03_Quantification$" = key_quantification, "04_Gene_Sets$" = key_gene_sets,
  "05_Set_Scores$" = key_set_scores, "Module_Building$" = key_module_building,
  "Module_Characterisation$" = key_characterisation, "Design_Contrasts$" = key_design,
  "Protein_Contrasts$" = key_protein, "Set_Contrasts$" = key_sets,
  "Module_Contrasts$" = key_modules, "Module_Preservation$" = key_preservation,
  "_groups/.*06_Classification/" = key_classification,
  "_groups/(.*/)?07_Trait_Associations$" = key_traits, "Clustering$" = key_clustering,
  "Cluster_Phenotype$" = key_cluster_phenotype,
  "01_Trait_Associations/" = key_trait_associations, "Trait_Pathways$" = key_trait_pathways,
  "Profile_Clusters$" = key_profiles
)

# Which shared function in R/ a step calls, if any.
code_of <- c(
  "Design_Contrasts$" = "R/design_contrasts.R", "Protein_Contrasts$" = "R/protein_contrasts.R",
  "Set_Contrasts$" = "R/set_contrasts.R", "Module_Contrasts$" = "R/module_contrasts.R",
  "06_Classification/" = "R/classify.R",
  "_groups/(.*/)?07_Trait_Associations$" = "R/group_trait_associations.R",
  "01_Trait_Associations/" = "R/associate.R"
)
first_match <- function(stage, table) {
  hit <- names(table)[map_lgl(names(table), \(pattern) grepl(pattern, stage))]
  if (length(hit)) table[[hit[1]]]
}

# The cluster groupings' mirrored steps and the HR/LR step whose method they share.
template_of <- c(
  "02_Design_Contrasts" = "01_Design_Contrasts", "03_Protein_Contrasts" = "02_Protein_Contrasts",
  "04_Set_Contrasts" = "03_Set_Contrasts", "05_Module_Contrasts" = "04_Module_Contrasts",
  "06_Classification/01_Proteins" = "06_Classification/01_Proteins",
  "06_Classification/02_Sets" = "06_Classification/02_Sets",
  "06_Classification/03_Modules" = "06_Classification/03_Modules",
  "07_Trait_Associations" = "07_Trait_Associations"
)
groups_of <- c(
  Muscle_size = "Hi and Lo, the muscle-size clusters",
  Muscle_strength = "Hi and Lo, the strength clusters",
  Proteome = "C1 and C2, the proteome clusters at each input"
)
level_words <- list(
  "01_Proteins" = c(noun = "protein", unit = "protein", collection = ""),
  "02_Sets" = c(noun = "gene-set score", unit = "set", collection = ", and within collection"),
  "03_Modules" = c(noun = "module eigengene", unit = "module", collection = "")
)

text_for <- function(stage) {
  key <- if (grepl("/06_Classification/", stage)) {
    "02_HRvLR_groups/06_Classification"
  } else if (grepl("^04_Continuous_no_group/01_Trait_Associations/", stage)) {
    "04_Continuous_no_group/01_Trait_Associations"
  } else {
    stage
  }
  text <- readme_text[[key]]
  words <- level_words[[basename(stage)]]
  if (is.null(words)) {
    return(text)
  }
  fill <- function(x) {
    for (w in names(words)) x <- gsub(paste0("{", w, "}"), words[[w]], x, fixed = TRUE)
    x
  }
  map(text, fill)
}

facts_of <- function(stage) {
  per_input <- grepl("^03_Cluster_groups/Proteome/0[2-7]", stage)
  first <- if (per_input) file.path(stage, "T1") else stage
  key_fun <- first_match(stage, key_of)
  key <- if (is.null(key_fun)) {
    character()
  } else if (per_input) {
    unlist(map(inputs_used, \(i) {
      c(sprintf("- `%s`", i), paste0("  ", key_fun(file.path(stage, i))))
    }))
  } else {
    key_fun(stage)
  }
  books <- list.files(here(first, "c_data"), "[.]xlsx$", full.names = TRUE)
  reads <- unique(unlist(map(books, \(b) rd(b, "input_manifest")$path)))
  if (per_input) reads <- gsub("/T1/", "/<input>/", gsub("/T1[.]rds$", "/<input>.rds", reads))
  sheets <- map(books, \(b) {
    rd(b, "read_me") |>
      filter(!sheet %in% c("columns", "input_manifest", "package_versions")) |>
      mutate(book = basename(b))
  }) |>
    list_rbind()
  listing <- function(dir) {
    files <- unlist(map(c("c_data", "b_reports"), \(d) {
      x <- list.files(here(dir, d))
      x <- x[!dir.exists(here(dir, d, x))]
      if (length(x)) file.path(d, x)
    }))
    files[!grepl("[.]html$", files)]
  }
  writes <- listing(stage)
  if (per_input) writes <- c(file.path("<input>", listing(first)), writes)
  list(
    key = key, reads = reads, writes = writes, sheets = sheets,
    script = file.path(stage, "a_script", list.files(here(stage, "a_script"), "[.]R$"))
  )
}

ticks <- function(x) paste0("`", x, "`", collapse = ", ")
files_table <- function(facts, code, note = "") {
  c(
    "| | |", "|---|---|",
    sprintf("| Reads | %s |", ticks(facts$reads)),
    sprintf("| Writes | %s |", ticks(facts$writes)),
    sprintf("| Run | `Rscript %s` |", facts$script),
    if (!is.null(code)) sprintf("| Code | `%s`%s |", code, note),
    ""
  )
}
sheets_table <- function(sheets) {
  if (!nrow(sheets)) {
    return(character())
  }
  books <- unique(sheets$book)
  body <- if (length(books) > 1) {
    c(
      "| Workbook | Sheet | Holds |", "|---|---|---|",
      sprintf("| `%s` | `%s` | %s |", sheets$book, sheets$sheet, sheets$holds)
    )
  } else {
    c(
      sprintf(
        paste(
          "`%s` opens with `read_me` and `columns` (every column defined) and closes with",
          "`input_manifest` and `package_versions`."
        ),
        books
      ),
      "", "| Sheet | Holds |", "|---|---|",
      sprintf("| `%s` | %s |", sheets$sheet, gsub("|", "\\|", sheets$holds, fixed = TRUE))
    )
  }
  c("## Sheets", "", body, "")
}
key_section <- function(key) if (length(key)) c("## Key result", "", key, "")

full_readme <- function(stage, facts) {
  text <- text_for(stage)
  code <- first_match(stage, code_of)
  note <- if (startsWith(stage, "02_HRvLR_groups") && !is.null(code) && code != "R/associate.R") {
    ", shared with the cluster groupings"
  } else {
    ""
  }
  c(
    paste("#", basename(stage)), "", text$question, "",
    files_table(facts, code, note),
    key_section(facts$key),
    "## Method", "", paste("-", text$method), "",
    if (length(text$care)) c("## Read with care", "", paste("-", text$care), ""),
    sheets_table(facts$sheets)
  )
}

mirror_readme <- function(stage, facts) {
  grouping <- strsplit(stage, "/")[[1]][2]
  template <- file.path("02_HRvLR_groups", template_of[[sub("^[^/]+/[^/]+/", "", stage)]])
  code <- first_match(template, code_of)
  up <- strrep("../", length(strsplit(stage, "/")[[1]]))
  how <- sprintf(
    paste(
      "The method is the HR/LR step's ([`%s`](%s%s/README.md)), run through the same function",
      "in `%s`."
    ),
    template, up, template, code
  )
  if (grouping == "Proteome") {
    how <- paste(
      how, "It runs once per clustering input, each into its own `<input>/` folder;",
      "`stack_inputs()` (`../R/inputs.R`) then joins the five workbooks into one, adding `input`",
      "and `validity` (circular, caution or clean) to every row."
    )
  }
  c(
    paste("#", basename(stage)), "", text_for(template)$question, "",
    sprintf("Groups: %s.", groups_of[[grouping]]), "", how, "",
    files_table(facts, code),
    key_section(facts$key),
    sheets_table(facts$sheets)
  )
}

for (stage in stages[!grepl("/Figures/", stages)]) {
  facts <- facts_of(stage)
  mirrored <- startsWith(stage, "03_Cluster_groups/") &&
    sub("^[^/]+/[^/]+/", "", stage) %in% names(template_of)
  lines <- if (mirrored) mirror_readme(stage, facts) else full_readme(stage, facts)
  while (lines[length(lines)] == "") lines <- lines[-length(lines)]
  writeLines(lines, here(stage, "README.md"))
}
