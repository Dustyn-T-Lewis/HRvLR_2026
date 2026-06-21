# F08 Acute Enrichment — Figure-specific style
# Enrichment analysis for Acute contrasts: Acute_HR, Acute_LR, Acute_Interaction
# Mirrors F09 concordance style but adds Interaction category
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/concordance_utils.R")

# Significance colors for Acute enrichment (including Interaction)
SIG_COLORS_F08 <- c(
  "Sig Both"    = "#2E7D32",
  "Sig HR only" = "#67A9CF",
  "Sig LR only" = "#EF8A62",
  "Interaction" = "#FF8F00",
  "NS"          = "grey70"
)
SIG_LABEL_FILL_F08 <- c(
  "Sig Both"    = scales::alpha("#2E7D32", 0.80),
  "Sig HR only" = scales::alpha("#67A9CF", 0.80),
  "Sig LR only" = scales::alpha("#EF8A62", 0.80),
  "Interaction" = scales::alpha("#FF8F00", 0.80),
  "NS"          = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F08 <- c(
  "Sig Both"    = "white",
  "Sig HR only" = "white",
  "Sig LR only" = "white",
  "Interaction" = "white",
  "NS"          = "white"
)

ORA_QUAD_COLORS_F08 <- c(
  "Concordant Up"      = "#D6604D",
  "Concordant Down"    = "#4393C3",
  "Discordant HR Up"   = "#7B5EA7",
  "Discordant LR Up"   = "#E6AB02"
)

# Classify proteins including Acute_Interaction
classify_proteins_f08 <- function(pi_hr, pi_lr, pi_int, threshold = 0.05) {
  dplyr::case_when(
    pi_hr < threshold & pi_lr < threshold ~ "Sig Both",
    pi_int < threshold                    ~ "Interaction",
    pi_hr < threshold                     ~ "Sig HR only",
    pi_lr < threshold                     ~ "Sig LR only",
    TRUE                                  ~ "NS"
  )
}
