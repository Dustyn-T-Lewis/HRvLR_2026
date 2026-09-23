# Pathway-level classification and phenotype association on singscore scores.
#
# The screens run once per collection so BH stays inside a collection, matching
# the set tests: a GO-Slim term and a GO:BP leaf are not the same size of
# question.

pacman::p_load(here, dplyr, purrr)

source(here("functions", "classify.R"))
source(here("functions", "shared_utils.R"))

OUT_DIR <- here("03_Pathways", "04_Classify_Associate", "c_data")

ss <- readRDS(here("03_Pathways", "03_Set_Scores", "c_data", "set_scores.rds"))
by_collection <- split(ss$catalog$set, ss$catalog$collection)

per_collection <- imap(by_collection, function(sets, collection) {
  screen_level(ss$scores[sets, , drop = FALSE], ss$meta) |>
    map(\(d) mutate(d, collection = collection, .before = 1))
})
screens <- set_names(names(per_collection[[1]])) |>
  map(\(part) list_rbind(map(per_collection, part)))
themes <- setNames(ss$catalog$theme, ss$catalog$set)
screens[c("classify", "associate")] <- lapply(
  screens[c("classify", "associate")],
  \(d) mutate(d, theme = themes[.data$feature], .after = "feature")
)

clear_dir(OUT_DIR)
saveRDS(screens, file.path(OUT_DIR, "pathway_screens.rds"))
openxlsx::write.xlsx(screens, file.path(OUT_DIR, "pathway_screens.xlsx"))

print(as.data.frame(screens$chance_classify), row.names = FALSE, digits = 3)
