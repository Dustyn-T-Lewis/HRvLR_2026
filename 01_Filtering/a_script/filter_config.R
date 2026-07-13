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
# proteins. There is a clean gap: every genuine plasma protein sits below myonuclei 1.
# Chosen on HPA annotation alone - no contrast, phenotype or group label was consulted.

filter_cfg <- list(
  hpa_file        = "HPA_annotations_full.tsv",
  miss_min_reps   = 5, # min detected samples in a Group_Time cell
  miss_min_groups = 1, # min Group_Time cells clearing miss_min_reps
  outlier_k       = 3, # methods that must agree for consensus (>=3/4)
  mahal_p         = 0.01, # PCA Mahalanobis chi-square tail cutoff
  mad_k           = 3, # Hampel constant for MAD-based flags
  ery_cut         = 5000, # erythrocyte nCPM at/above = red-cell protein
  myo_cut         = 20, # myonuclei nCPM at/above = candidate for muscle rescue
  blood_max       = 1e9 # rescue only below this blood conc (pg/L); above = true plasma
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
