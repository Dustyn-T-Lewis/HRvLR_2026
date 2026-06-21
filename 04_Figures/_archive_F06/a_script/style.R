# F06 Chronic vs Acute — Figure-specific style
# Compares Training (chronic adaptation) vs Acute exercise within each group
# Two parallel figures: F06_HR (Training_HR vs Acute_HR), F06_LR (Training_LR vs Acute_LR)
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/concordance_utils.R")

# Significance colors — HR panels
SIG_COLORS_F06_HR <- c(
  "Sig Both"          = "#2E7D32",
  "Sig Training only" = "#2166AC",
  "Sig Acute only"    = "#67A9CF",
  "NS"                = "grey70"
)
SIG_LABEL_FILL_F06_HR <- c(
  "Sig Both"          = scales::alpha("#2E7D32", 0.80),
  "Sig Training only" = scales::alpha("#2166AC", 0.80),
  "Sig Acute only"    = scales::alpha("#67A9CF", 0.80),
  "NS"                = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F06_HR <- c(
  "Sig Both"          = "white",
  "Sig Training only" = "white",
  "Sig Acute only"    = "white",
  "NS"                = "white"
)

# Significance colors — LR panels
SIG_COLORS_F06_LR <- c(
  "Sig Both"          = "#2E7D32",
  "Sig Training only" = "#B2182B",
  "Sig Acute only"    = "#EF8A62",
  "NS"                = "grey70"
)
SIG_LABEL_FILL_F06_LR <- c(
  "Sig Both"          = scales::alpha("#2E7D32", 0.80),
  "Sig Training only" = scales::alpha("#B2182B", 0.80),
  "Sig Acute only"    = scales::alpha("#EF8A62", 0.80),
  "NS"                = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F06_LR <- c(
  "Sig Both"          = "white",
  "Sig Training only" = "white",
  "Sig Acute only"    = "white",
  "NS"                = "white"
)

# Concordance quadrant colors (same for both groups)
ORA_QUAD_COLORS_F06 <- c(
  "Concordant Up"   = "#D6604D",
  "Concordant Down" = "#4393C3",
  "Discordant Training Up" = "#7B5EA7",
  "Discordant Acute Up"    = "#E6AB02"
)
