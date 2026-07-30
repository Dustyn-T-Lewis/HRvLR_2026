# Filtering cutoffs and the curated contaminant list, single-sourced by the run script
# and the filtering tutorial so the two cannot drift.
#
# HPA source: 00_input/HPA_annotations_full.tsv, retrieved 2026-07-13 from
#   proteinatlas.org/api/search_download.php?search=&format=tsv&compress=no
#   &columns=g,up,pc,secl,blconcms,sc_RNA_Myonuclei,sc_RNA_Erythrocytes
# All 20162 genes, no tissue filter. It is an annotation lookup, not a protein universe:
# our universe is the 2400 proteins the MS measured. The previous 12893-row export was a
# muscle-expressed query (94% of its rows have myonuclei RNA > 0, vs 51% of the rows it
# omitted), so gating on "present in HPA" silently discarded any measured protein whose
# transcript the single-cell atlas did not see in a myonucleus - a transcriptomic prior on
# a proteomic dataset. It took CKMT1A and ATP1A3 (both myonuclei nCPM 0) with it.
#
# myo_cut is 20, not 50. HPA tags any plasma-detectable protein "Secreted to blood", which
# is why CK and myoglobin are clinical damage markers; the rescue exists to undo that. At 50
# it failed to fire for GPI (myonuclei 43.1), ANXA2 (40.2) and PPIA (21.1), all real muscle
# proteins. Chosen on HPA annotation alone - no contrast, phenotype or group label was
# consulted, so the choice cannot leak into the HR/LR result.
#
# There is NO clean gap on myonuclei, and an earlier version of this comment claimed there was.
# Removed-as-plasma proteins span myonuclei 0-236.8 and rescued ones span 20.6-362.3; they
# overlap across the whole rescue range, and 8 removed proteins sit above myo_cut (BTD 236.8,
# A2M 105.1, C1S 89.0, PZP 74.7, ITIH2 44.0, SERPINF1 28.3, GSN 26.0, ITIH4 21.4). blood_max
# is what actually separates them: all 8 measure 2.3e9-4.2e11 pg/L in plasma, while genuinely
# intracellular proteins measure 1e6-1e8. The measured concentration does the work; the
# transcript floor only decides who gets a second hearing.
#
# The floor is therefore imperfect, and it is known to over-delete: C1QBP (mitochondrial,
# blood 1.5e6), HMGB2 (nuclear, 8.3e5), PPIB, LGALS1 and CTSB (an exercise myokine) all clear
# blood_max and fail only on myonuclei. 03_DEP/a_non_imputed/a_script/02_blood_filter_sensitivity.R
# refits the nine contrasts with all 136 blood-tagged proteins readmitted: zero BH hits in every
# HR-vs-LR and interaction contrast, so the null does not depend on this rule. The same run shows
# why the rule is kept - readmitting them floods Acute_LR from 19 to 97 hits, 57 of them the
# readmitted blood proteins, which is the 2x-bloodier T3 biopsy confound arriving on cue.

# blood_cor is a reported diagnostic, not a gate. The index is the per-sample
# mean log2 intensity of the haemoglobin anchor, measured before removal, and
# each protein's Spearman correlation with it says how far that protein tracks
# carryover in THIS dataset. It travels in protein_calls.csv and F04 flags rows
# by it.
#
# It stopped gating removal on 2026-07-30. The cut was 0.45, justified as a
# rounding above a permutation null recorded as 0.43; recomputing that null gave
# 0.59, and stratifying it showed the quantile depends on how many samples a
# protein was seen in - 0.71 at 15-25 observations against 0.47 at 41-48. One
# flat cut cannot serve both ends of that range. Removal moved to the curated
# blood list below, whose membership is protein identity.
#
# rbc_file still corroborates: it is what makes the enucleate red-cell membrane
# skeleton visible at all, since band 3, spectrin, protein 4.2 and CA1 are
# RNA-poor and escape ery_cut. It never removes on its own, because a deep RBC
# proteome shares roughly 70% of any muscle proteome.
filter_cfg <- list(
  hpa_file        = "HPA_annotations_full.tsv",
  rbc_file        = "RBC_proteome_reference.tsv", # five-source mature-RBC proteome, built in stage 00
  miss_min_reps   = 5, # min detected samples in a Group_Time cell
  miss_min_groups = 1, # min Group_Time cells clearing miss_min_reps
  outlier_k       = 3, # methods that must agree for consensus (>=3/4)
  mahal_p         = 0.01, # PCA Mahalanobis chi-square tail cutoff
  mad_k           = 3, # Hampel constant for MAD-based flags
  ery_cut         = 5000, # erythrocyte nCPM at/above = red-cell protein
  myo_cut         = 20, # myonuclei nCPM at/above = candidate for muscle rescue
  blood_max       = 1e9, # rescue only below this blood conc (pg/L); above = true plasma
  blood_anchor    = c("HBB", "HBA1", "HBD", "HBG1", "HBG2") # haemoglobins define the index
)

# Handling contaminants no annotation rule can catch. Matched by ACCESSION only: cRAP
# collides with real muscle proteins on gene symbol (ALDOA vs rabbit ALDOA_RABIT; GLUD1 vs
# bovine P00366, one digit from human P00367) and on description (a /trypsin|casein|albumin/
# regex eats parvalbumin and the casein kinases). Only the cRAP dust/contact section is used
# - the Sigma UPS1 section is a spike-in quantitation standard we never spiked, and it
# contains myoglobin and creatine kinase M, ranks 1 and 2 of this proteome.
contaminants <- tibble::tribble(
  ~uniprot_id, ~gene,    ~class,           ~reason,
  "P04264",    "KRT1",   "keratin",        "cRAP dust/contact; cornified epidermal",
  "P35908",    "KRT2",   "keratin",        "cRAP dust/contact; epidermal",
  "P13645",    "KRT10",  "keratin",        "cRAP dust/contact; type-I partner of KRT1",
  "P35527",    "KRT9",   "keratin",        "cRAP dust/contact; palmoplantar",
  "P08779",    "KRT16",  "keratin",        "cRAP dust/contact",
  "P04259",    "KRT6B",  "keratin",        "cRAP dust/contact",
  "P02533",    "KRT14",  "keratin",        "cRAP dust/contact; basal epidermal",
  "P13647",    "KRT5",   "keratin",        "cRAP dust/contact; basal epidermal",
  "P02538",    "KRT6A",  "keratin",        "cRAP dust/contact",
  "Q04695",    "KRT17",  "keratin",        "cRAP dust/contact",
  "Q7Z794",    "KRT77",  "keratin",        "cRAP dust/contact",
  "Q8N1N4",    "KRT78",  "keratin",        "cRAP dust/contact",
  "P05787",    "KRT8",   "keratin",        "simple-epithelial, but Spearman +0.78 with the KRT1/2/10 load across samples",
  "P15924",    "DSP",    "keratin",        "desmoplakin; skin desmosome",
  "Q02413",    "DSG1",   "keratin",        "desmoglein-1; skin desmosome",
  "P81605",    "DCD",    "keratin",        "dermcidin; sweat",
  "P31151",    "S100A7", "keratin",        "psoriasin; skin",
  "P01040",    "CSTA",   "keratin",        "cystatin-A; cornified envelope",
  "P68871",    "HBB",    "globin",         "red-cell carryover",
  "P69905",    "HBA1",   "globin",         "red-cell carryover",
  "P02042",    "HBD",    "globin",         "red-cell carryover",
  "P69891",    "HBG1",   "globin",         "red-cell carryover",
  "P69892",    "HBG2",   "globin",         "orphan gamma-globin: HbF is a2g2 and no alpha or beta chain survives. Shares tryptic peptides with HBB/HBD; razor-assigned onto the last entry standing",
  "P09105",    "HBQ1",   "globin",         "red-cell carryover",
  "Q9NZD4",    "AHSP",   "globin",         "erythroid-specific haemoglobin chaperone",
  "P01834",    "IGKC",   "immunoglobulin", "absent from HPA; constant region",
  "P01860",    "IGHG3",  "immunoglobulin", "absent from HPA; constant region",
  "P01880",    "IGHD",   "immunoglobulin", "absent from HPA; constant region",
  "P02746",    "C1QB",   "complement",     "absent from HPA",
  "P00736",    "C1R",    "complement",     "absent from HPA"
)

# The blood list is curated by protein identity, not by covariation. Every member
# is an immunoglobulin, a complement component, a secreted plasma protein, a
# mature red-cell protein or a leukocyte protein, and carries that reason in the
# file. It replaces the blood_cor > 0.45 removal rule, whose 0.43 justification
# did not reproduce and whose null turned out to depend on how many samples a
# protein was seen in (0.71 at 15-25 observations against 0.47 at 41-48), so one
# flat cut over-removed sparse proteins and under-removed dense ones.
#
# What that rule caught and this list does not: 16 nuclear and cytoskeletal
# proteins with real myonuclei expression, RBMX at 358 nCPM and SYNE2 at 303
# among them, where covariation with haemoglobin was the only evidence. They are
# back in the matrix. See 00_input/PROVENANCE.md.
contaminants <- dplyr::bind_rows(
  contaminants,
  readr::read_csv(
    here::here("00_input", "blood_contaminants.csv"),
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    dplyr::select("uniprot_id", "gene", "class", "reason")
)

# Reported per sample as % of summed intensity, measured BEFORE removal. These are
# sample-quality readouts, not filters: an imbalance across HR/LR or timepoint is a protocol
# finding, not a nuisance to delete.
qc_panels <- list(
  keratin   = c("KRT1", "KRT2", "KRT10", "KRT8", "KRT9", "KRT16", "KRT14", "KRT5", "KRT6A", "KRT6B", "DCD"),
  blood     = c("HBB", "HBA1", "HBD", "HBG1", "HBG2", "HBQ1", "AHSP", "ALB", "CA1", "SLC4A1", "SPTA1"),
  leukocyte = c("MPO", "CTSG", "ELANE", "PRTN3", "AZU1", "S100A8", "S100A9", "LYZ", "LTF"),
  adipose   = c("FABP4", "PLIN1", "PLIN4", "ADIPOQ", "ADIRF"),
  myofibre  = c("MB", "CKM", "ACTA1", "MYH1", "MYH2", "MYH7", "ALDOA", "CASQ1", "PYGM", "ATP2A1", "DES")
)
