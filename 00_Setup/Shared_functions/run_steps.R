# Render one grouping's step scripts: <phase>/code/<step>.R becomes
# <phase>/<grouping>/b_reports/<step>.pdf, with params$grouping set. Run from the project root,
# so Quarto's R session starts with the project's renv library.
run_steps <- function(phase, grouping, steps) {
  for (step in steps) {
    status <- system2("quarto", c(
      "render", file.path(phase, "code", paste0(step, ".R")),
      "-P", paste0("grouping:", grouping),
      "--output-dir", file.path("..", grouping, "b_reports")
    ))
    if (status != 0) stop(phase, "/code/", step, ".R failed for ", grouping)
  }
}
