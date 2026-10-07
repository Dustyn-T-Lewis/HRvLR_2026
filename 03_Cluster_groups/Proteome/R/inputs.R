# The five proteome clusterings, the time points each was built from, and the validity of every
# comparison against them. A comparison that uses exactly the clustering's time points is
# circular (the clustering made the difference); one that shares a time point is caution
# (regression to the mean can make a difference); the rest are clean.

inputs_used <- c("T1", "T2", "T3", "training", "acute")
time_points <- list(
  T1 = "T1", T2 = "T2", T3 = "T3", training = c("T1", "T2"), acute = c("T2", "T3"),
  all = c("T1", "T2", "T3")
)

# The time points a contrast, task or window name uses; PERMANOVA's windows name a time point, or
# every time point ("all", "subject mean").
uses <- function(name) {
  dplyr::case_when(
    name %in% c("T1", "T2", "T3") ~ name,
    name %in% c("all", "subject mean") ~ "all",
    grepl("^Acute_C1vC2$|^post_bout_level$", name) ~ "T3",
    grepl("^Baseline|^baseline", name) ~ "T1",
    grepl("^Trained|^post_level$", name) ~ "T2",
    grepl("^Training|^training$", name) ~ "training",
    grepl("^Acute|^acute$", name) ~ "acute"
  )
}
validity <- function(input, name) {
  own <- time_points[[input]]
  purrr::map_chr(uses(name), \(used) {
    if (is.na(used)) {
      return(NA_character_)
    }
    points <- time_points[[used]]
    if (setequal(points, own)) "circular" else if (any(points %in% own)) "caution" else "clean"
  })
}

# One workbook, <content>_all_inputs.xlsx, from the five per-input workbooks of a stage: every
# result sheet stacked with an
# input column, and a validity column wherever a contrast, task or window names the comparison.
# Sheets are read as text and typed after stacking, since an empty sheet reads as logical. Sheets in
# `leave_out` stay in the per-input workbooks only (xlsx caps a sheet at 1,048,576 rows).
stack_inputs <- function(stage, leave_out = character()) {
  files <- purrr::map_chr(inputs_used, \(input) {
    list.files(file.path(stage, input, "c_data"), "\\.xlsx$", full.names = TRUE)
  })
  skip <- c("read_me", "input_manifest", "package_versions", "plot_index", leave_out)
  sheets <- setdiff(readxl::excel_sheets(files[1]), skip)
  stacked <- purrr::map(rlang::set_names(sheets), \(sheet) {
    purrr::map2(files, inputs_used, \(file, input) {
      readxl::read_excel(file, sheet, col_types = "text") |>
        dplyr::mutate(input = input, .before = 1)
    }) |>
      purrr::list_rbind() |>
      readr::type_convert(col_types = readr::cols())
  }) |>
    purrr::map(\(d) {
      key <- intersect(c("contrast", "task", "window"), names(d))[1]
      if (is.na(key)) {
        return(d)
      }
      dplyr::mutate(d, validity = validity(input[1], .data[[key]]), .by = input, .after = input)
    })
  dir.create(file.path(stage, "c_data"), showWarnings = FALSE)
  write_workbook(
    stacked,
    paste(
      "Sheet", names(stacked), "from every clustering input, stacked; validity is circular,",
      "caution or clean against that input."
    ),
    file.path(stage, "c_data", sub("[.]xlsx$", "_all_inputs.xlsx", basename(files[1]))),
    rlang::set_names(sub(paste0(here::here(), "/"), "", files), inputs_used), files,
    note = if (length(leave_out)) {
      paste("Not stacked:", paste(leave_out, collapse = ", "), "- see each input's own workbook.")
    }
  )
}
