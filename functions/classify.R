# Can one feature tell two sets of samples apart?
#
# Seven tasks mirror the contrasts. Four ask whether a feature separates a later
# timepoint from an earlier one inside one arm, pairing each subject with
# itself. Three ask whether it separates HR from LR, on the baseline level or on
# a subject's own change, one value per subject.
#
# The AUC is read as a direction as well as a size: above 0.5 means higher in
# the case group (the later timepoint, or HR). The p comes from the Wilcoxon
# test the AUC is a rescaling of, signed-rank when the samples pair. BH runs
# within a task, never across tasks or feature levels, because each task is its
# own question with its own null.

pacman::p_load(here, dplyr, tibble, purrr, pROC)

source(here("functions", "association.R"))

TASKS <- tribble(
  ~task,            ~compare, ~arm, ~control, ~case,
  "Training_HR",    "time",   "HR", "T1",     "T2",
  "Training_LR",    "time",   "LR", "T1",     "T2",
  "Acute_HR",       "time",   "HR", "T2",     "T3",
  "Acute_LR",       "time",   "LR", "T2",     "T3",
  "Baseline_HRvLR", "arm",    NA,   "LR",     "HR",
  "Training_HRvLR", "arm",    NA,   "LR",     "HR",
  "Acute_HRvLR",    "arm",    NA,   "LR",     "HR"
)

# Which subject_window() an arm task reads.
ARM_WINDOW <- c(
  Baseline_HRvLR = "T1", Training_HRvLR = "training", Acute_HRvLR = "acute"
)

MIN_PER_GROUP <- 3L

auc_one <- function(control, case, paired) {
  keep <- if (paired) !is.na(control) & !is.na(case) else TRUE
  control <- control[keep & !is.na(control)]
  case <- case[keep & !is.na(case)]
  if (min(length(control), length(case)) < MIN_PER_GROUP) {
    return(c(auc = NA_real_, p = NA_real_))
  }
  roc <- pROC::roc(
    controls = control, cases = case, direction = "<", quiet = TRUE
  )
  test <- stats::wilcox.test(case, control, paired = paired, exact = FALSE)
  c(auc = as.numeric(pROC::auc(roc)), p = test$p.value)
}

# Two feature-by-subject matrices, columns in matching subject order.
task_matrices <- function(mat, task, meta) {
  arm_of <- distinct(meta, .data$subject, .data$arm)
  arm_of <- setNames(arm_of$arm, arm_of$subject)
  if (task$compare == "time") {
    control <- subject_window(mat, task$control, meta)
    case <- subject_window(mat, task$case, meta)
    both <- intersect(colnames(control), colnames(case))
    both <- both[arm_of[both] == task$arm]
    return(list(control = control[, both], case = case[, both]))
  }
  x <- subject_window(mat, ARM_WINDOW[[task$task]], meta)
  list(
    control = x[, arm_of[colnames(x)] == task$control, drop = FALSE],
    case = x[, arm_of[colnames(x)] == task$case, drop = FALSE]
  )
}

classify_features <- function(mat, task, meta = sample_metadata()) {
  m <- task_matrices(mat, task, meta)
  paired <- task$compare == "time"
  stats <- vapply(
    seq_len(nrow(mat)),
    function(i) auc_one(m$control[i, ], m$case[i, ], paired),
    numeric(2)
  )
  tibble(
    task = task$task,
    feature = rownames(mat),
    n_control = ncol(m$control),
    n_case = ncol(m$case),
    auc = stats["auc", ],
    p = stats["p", ]
  ) |>
    mutate(bh = stats::p.adjust(.data$p, "BH"))
}

classify_all <- function(mat, meta = sample_metadata(), tasks = TASKS) {
  map(seq_len(nrow(tasks)), \(i) classify_features(mat, tasks[i, ], meta)) |>
    list_rbind()
}

# Both screens and their chance tables for one feature level, the unit every
# stage's classify-and-associate step writes.
screen_level <- function(mat, meta = sample_metadata(),
                         pheno = phenotype_table()) {
  classified <- classify_all(mat, meta)
  associated <- associate_all(mat, meta, pheno)
  list(
    classify = classified,
    associate = associated,
    chance_classify = chance_table(classified, "task"),
    chance_associate = chance_table(associated, "window", "phenotype")
  )
}
