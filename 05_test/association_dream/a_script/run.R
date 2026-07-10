# Run the association_dream unit end to end from the repo root.
pacman::p_load(here)
a <- here("05_test", "association_dream", "a_script")
source(file.path(a, "setup.R"))
source(file.path(a, "01_dream_continuous.R"))
source(file.path(a, "02_dream_categorical.R"))
source(file.path(a, "03_figures.R"))
