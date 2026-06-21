# Concordance / Acute-change analysis utilities
# Shared across F4 (between-group) and F5 (within-group) figures.

# --- Significance classification for concordance scatter ---
classify_concordance <- function(pi_x, pi_y, label_x, label_y, threshold = 0.05) {
  dplyr::case_when(
    pi_x < threshold & pi_y < threshold ~ "Sig Both",
    pi_x < threshold                    ~ paste("Sig", label_x, "only"),
    pi_y < threshold                    ~ paste("Sig", label_y, "only"),
    TRUE                                ~ "NS"
  )
}

# --- Concordance color palette generator ---
make_concordance_colors <- function(label_x, label_y) {
  c(
    "Sig Both"                          = "#2E7D32",
    setNames("#4393C3", paste("Sig", label_x, "only")),
    setNames("#D6604D", paste("Sig", label_y, "only")),
    "NS"                                = "grey70"
  )
}

# --- Quadrant classification ---
classify_quadrant <- function(lfc_x, lfc_y, label_x, label_y) {
  dplyr::case_when(
    lfc_x > 0 & lfc_y > 0 ~ "Concordant Up",
    lfc_x < 0 & lfc_y < 0 ~ "Concordant Down",
    lfc_x > 0 & lfc_y < 0 ~ paste0("Discordant (", label_x, " Up / ", label_y, " Down)"),
    lfc_x < 0 & lfc_y > 0 ~ paste0("Discordant (", label_x, " Down / ", label_y, " Up)"),
    TRUE                   ~ NA_character_
  )
}

QUAD_COLORS <- c(
  "Concordant Up"   = "#E57373",
  "Concordant Down" = "#64B5F6"
)

# --- Concordance statistics ---
concordance_stats <- function(lfc_x, lfc_y, pi_x = NULL, pi_y = NULL, threshold = 0.05) {
  valid <- !is.na(lfc_x) & !is.na(lfc_y)
  x <- lfc_x[valid]; y <- lfc_y[valid]
  n <- length(x)

  # All proteins
  r_all  <- cor(x, y, method = "pearson")
  rho_all <- cor(x, y, method = "spearman")
  sign_conc_all <- mean(sign(x) == sign(y)) * 100
  r_ci   <- fisher_z_ci(r_all, n)
  rho_ci <- fisher_z_ci(rho_all, n)

  stats <- list(
    n_all = n,
    pearson_all = r_all, pearson_ci = r_ci,
    spearman_all = rho_all, spearman_ci = rho_ci,
    sign_concordance_all = sign_conc_all
  )

  # Sig-only (if pi-scores provided)
  if (!is.null(pi_x) && !is.null(pi_y)) {
    sig <- valid & (!is.na(pi_x) & pi_x < threshold | !is.na(pi_y) & pi_y < threshold)
    xs <- lfc_x[sig]; ys <- lfc_y[sig]; ns <- length(xs)
    if (ns >= 5) {
      stats$n_sig <- ns
      stats$pearson_sig  <- cor(xs, ys, method = "pearson")
      stats$spearman_sig <- cor(xs, ys, method = "spearman")
      stats$sign_concordance_sig <- mean(sign(xs) == sign(ys)) * 100
    }
  }
  stats
}

# RRHO2 grid computation removed — use RRHO2::RRHO2_initialize() directly with
# stepsize = 20 (matches Cahill et al. 2018 and YvO reference). Callers:
# F05/panel_B_rrho2.R, F06/panel_B_rrho2.R, F09/supp/c02_RRHO2.R.
