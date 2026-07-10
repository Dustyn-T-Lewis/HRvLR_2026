# Run the singscore_vs_gsva_validation unit end to end from the repo root.
pacman::p_load(here)
a <- here("05_test", "singscore_vs_gsva_validation", "a_script")
source(file.path(a, "setup.R"))
source(file.path(a, "01_agreement.R"))
