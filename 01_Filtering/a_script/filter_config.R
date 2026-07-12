# Filtering cutoffs, single-sourced by the run script and the filtering tutorial
# so the two cannot drift. Values are HPA-grounded heuristics, not published
# thresholds; the blood classes follow the CvH design (secreted-to-blood,
# immunoglobulin, erythrocyte) with a myonuclei rescue.
#
# ery_cut, myo_cut, blood_max read HPA Single Cell Type RNA (nCPM) and blood MS
# concentration (pg/L): Uhlen et al. 2015 Science 347:1260419 (tissue atlas) and
# Uhlen et al. 2019 Science 366:eaax9198 (blood atlas). nCPM floors are set
# conservatively to catch red-cell and plasma contaminants while sparing genuine
# muscle proteins; blood_max vetoes the rescue for true high-abundance plasma
# proteins.

filter_cfg <- list(
  miss_min_reps   = 5, # min detected samples per Group_Time cell (8, or 7 where S28/S29 partial)
  miss_min_groups = 1, # min Group_Time cells clearing miss_min_reps
  outlier_k       = 3, # methods that must agree for consensus (>=3/4)
  mahal_p         = 0.01, # PCA Mahalanobis chi-square tail cutoff
  mad_k           = 3, # Hampel constant for MAD-based flags
  ery_cut         = 5000, # erythrocyte nCPM at/above = red-cell (hemoglobin) protein
  myo_cut         = 50, # myonuclei nCPM at/above = candidate for muscle rescue
  blood_max       = 1e9 # rescue only below this blood conc (pg/L); above = true plasma
)
