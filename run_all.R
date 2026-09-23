# Rebuild the pipeline from the filtered input to the three packets.
#
# Each step runs in its own R session, as YvO's runner does, so no object
# leaks from one stage into the next: a script that only works after another
# has left something in the workspace fails here instead of in a clean clone.
# 00_input's builders and the V1 equivalence check are left out; the first
# rewrites tracked inputs and the second needs the sibling V1 tree.
#
#   Rscript run_all.R

STEPS <- c(
  "01_Preprocess/01_Filtering/a_script/01_run_filtering.R",
  "01_Preprocess/01_Filtering/a_script/02_filtering_figure.R",
  "01_Preprocess/02_Normalization/a_script/01_run_normalization.R",
  "01_Preprocess/03_Imputation/a_script/01_impute_missforest.R",
  "02_Proteins/01_Differential/a_script/01_run_dep.R",
  "02_Proteins/01_Differential/a_script/03_fry_concordance.R",
  "02_Proteins/02_Classify_Associate/a_script/01_classify_associate.R",
  "02_Proteins/03_Packet/a_script/01_packet.R",
  "03_Pathways/01_Gene_Sets/a_script/01_build_gene_sets.R",
  "03_Pathways/02_Set_Tests/a_script/01_set_tests.R",
  "03_Pathways/03_Set_Scores/a_script/01_set_scores.R",
  "03_Pathways/04_Classify_Associate/a_script/01_classify_associate.R",
  "03_Pathways/05_Packet/a_script/01_packet.R",
  "04_Networks/01_Modules/a_script/01_build_modules.R",
  "04_Networks/02_Characterise/a_script/01_characterise_modules.R",
  "04_Networks/03_Classify_Associate/a_script/01_classify_associate.R",
  "04_Networks/04_Packet/a_script/01_packet.R"
)

log_dir <- here::here(".runlogs", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(log_dir, recursive = TRUE)

for (step in STEPS) {
  log <- file.path(log_dir, paste0(gsub("/", "__", step), ".log"))
  started <- Sys.time()
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--no-save", "--no-restore", shQuote(here::here(step))),
    stdout = log, stderr = log
  )
  message(sprintf(
    "%-70s %6.1fs", step, difftime(Sys.time(), started, units = "secs")
  ))
  if (status != 0) stop("failed: ", step, "\nlog: ", log, call. = FALSE)
}
message("done; logs in ", log_dir)
