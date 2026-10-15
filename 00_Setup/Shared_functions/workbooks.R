# Workbook writing for every step: each workbook opens with read_me and columns (every column
# described from workbook_columns.csv) and closes with the inputs' md5 and the loaded packages.

# What every column in every workbook holds, from 00_Setup/Shared_functions/workbook_columns.csv:
# a name (or a pattern starting with ^) and, where one name means different things, the sheets it
# applies to. The columns sheet lists each sheet's columns with their meaning; a column the
# dictionary does not describe stops the run, so no workbook ships an unexplained column.
column_dictionary <- readr::read_csv(
  here::here("00_Setup", "Shared_functions", "workbook_columns.csv"),
  col_types = "ccc", na = ""
)
describe_columns <- function(sheets) {
  rows <- purrr::imap(sheets, \(d, sheet) tibble(sheet = sheet, column = names(d))) |>
    purrr::list_rbind()
  meaning <- purrr::map2_chr(rows$sheet, rows$column, \(sheet, column) {
    fits <- column_dictionary |>
      filter(
        is.na(.data$sheet) | stringr::str_detect(.env$sheet, coalesce(.data$sheet, "^$")),
        if_else(startsWith(.data$column, "^"),
          stringr::str_detect(.env$column, .data$column), .data$column == .env$column
        )
      ) |>
      arrange(is.na(sheet))
    if (nrow(fits)) fits$meaning[1] else NA_character_
  })
  missing <- unique(rows$column[is.na(meaning)])
  if (length(missing)) {
    stop(
      "00_Setup/Shared_functions/workbook_columns.csv does not describe: ",
      paste(missing, collapse = ", ")
    )
  }
  mutate(rows, meaning = meaning)
}

# Every workbook opens on read_me and columns and ends with the inputs' md5 and the loaded packages;
# note, when given, is one more read_me row.
write_workbook <- function(sheets, holds, file, inputs, paths, note = NULL) {
  sheets <- c(sheets, list(
    input_manifest = tibble(
      input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
    ),
    package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
      as_tibble() |>
      select(package, version = loadedversion, source)
  ))
  holds <- c(holds, "Files read, with md5.", "Packages loaded at run time.")
  read_me <- tibble(sheet = c("columns", names(sheets)), holds = c(
    "Every column of every sheet below, with what it holds.", holds
  ))
  if (!is.null(note)) read_me <- add_row(read_me, sheet = "note", holds = note)
  columns <- describe_columns(sheets)
  writexl::write_xlsx(c(list(read_me = read_me, columns = columns), sheets), file)
}
