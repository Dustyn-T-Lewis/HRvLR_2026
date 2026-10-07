# The question, method and cautions of every step README, keyed by step folder. docs/readmes.R
# adds the files, the key result and the sheets from the saved outputs. In the classification and
# trait-association entries, {noun}, {unit} and {collection} are filled per level.
readme_text <- list(
  "01_Processing/01_Protein_Filtering" = list(
    question = "Which proteins and which runs are fit to analyse?",
    method = c(
      paste(
        "Contaminants go by identity, never by behaviour: the 30 accessions written into the",
        "script (keratins, globins, immunoglobulins, complement), the 95 in",
        "`00_Input/blood_contaminants.csv`, and proteins the Human Protein Atlas tags as blood",
        "unless their myonuclear RNA marks them as muscle."
      ),
      paste(
        "A run goes when three of four checks flag it: missingness, Mahalanobis distance on the",
        "first three principal components, median intensity (Hampel, 3 MAD) and median",
        "correlation with other runs."
      ),
      paste(
        "A protein stays when it is detected in at least 5 samples of at least one of the six",
        "group-by-time cells (`proteoDA::filter_proteins_by_group()`)."
      ),
      paste(
        "The blood index is the mean log2 intensity of HBB and HBA1 minus the sample's median",
        "log2 intensity (`blood_rel`). Later steps use it as a covariate or regress it out."
      )
    ),
    care = c(
      paste(
        "The detection filter is the one step in 01_Processing that reads the HR/LR label. It",
        "counts detections, never abundance."
      ),
      paste(
        "`blood_cor` and `in_rbc` are reported and decide nothing; `blood_cor_null` shows why",
        "no single correlation cut fits every protein."
      )
    )
  ),
  "01_Processing/02_Normalisation" = list(
    question = paste(
      "How are the 45 samples made comparable before any protein is compared across them?"
    ),
    method = c(
      paste(
        "log2, then cyclic loess (`proteoDA::normalize_data(norm_method = \"cycloess\")`).",
        "Missing values stay missing; `03_Quantification` models them."
      ),
      paste(
        "proteoDA's standard reports compare every normalisation method (`norm_comparison.pdf`)",
        "and show sample QC before (`qc_pre.pdf`) and after (`qc_post.pdf`)."
      ),
      paste(
        "`eta_squared` is the share of each protein's variance the six group-by-time cells",
        "explain (`limma::lmFit(~ group)`); descriptive, read by nothing downstream."
      )
    ),
    care = c(
      paste(
        "The PCA fills gaps with the sample median for display only; the saved matrix keeps",
        "them missing."
      )
    )
  ),
  "01_Processing/03_Quantification" = list(
    question = paste(
      "How does every protein get a value in every sample without treating estimates as",
      "measurements?"
    ),
    method = c(
      paste(
        "`limpa::dpc()` fits one detection-probability curve: the chance a value is detected,",
        "as a logistic function of its log2 intensity."
      ),
      paste(
        "`limpa::dpcQuantByRow()` gives each protein a value and a standard error in every",
        "sample. Observed values keep a small standard error; estimated ones get a wide one, so",
        "later fits weight them less."
      ),
      paste(
        "`DAList_corrected.rds` removes each protein's blood slope (`limma::removeBatchEffect`,",
        "design `~ 0 + timepoint`) for the methods that cannot take a covariate."
      )
    ),
    care = c(
      paste(
        "`prop_observed` in `proteins` says how much of a row was measured; a low share rests",
        "on the curve more than on data."
      )
    )
  ),
  "01_Processing/04_Gene_Sets" = list(
    question = paste(
      "Which groups of proteins count as one biological process, decided once for every",
      "analysis?"
    ),
    method = c(
      paste(
        "MSigDB 2026.1.Hs (Hallmark, Reactome, KEGG Legacy, GO:BP) from the snapshot in",
        "`00_Input`, plus GO Slim terms built from the GO:BP hierarchy (`GSEABase`, `GO.db`)."
      ),
      paste(
        "One protein per gene symbol: the one observed in the most samples. Proteins with no",
        "symbol or several are left out."
      ),
      paste(
        "A set is tested when it holds 15 to 500 genes in the source and at least 15 measured",
        "here (GO Slim: 15 to 500 measured)."
      )
    ),
    care = c(
      "Overlapping sets are kept, so one signal can appear in several neighbouring sets."
    )
  ),
  "01_Processing/05_Set_Scores" = list(
    question = "How strongly is each gene set expressed in each biopsy?",
    method = c(
      paste(
        "`singscore::rankGenes()` ranks the proteins within each sample;",
        "`singscore::multiScore()` gives one score per qualifying set and sample. Every set is",
        "scored as an up-set."
      ),
      paste(
        "`structure_check`: PCA of the samples on their scores and the share of each component",
        "that subject identity explains."
      )
    ),
    care = c(
      paste(
        "Scores follow the person more than the time point (`subject_share`), which is why the",
        "analyses block on subject and use within-person change."
      )
    )
  ),
  "01_Processing/06_Module_Building" = list(
    question = paste(
      "Which proteins rise and fall together within a person, and how does each such group move",
      "in each biopsy?"
    ),
    method = c(
      paste(
        "WGCNA (`blockwiseModules()`): signed network, bicor, `deepSplit = 2`, `minModuleSize =",
        "30`, `mergeCutHeight = 0.15`, `randomSeed = 54321`; soft power from",
        "`pickSoftThreshold()` (R² 0.85)."
      ),
      paste(
        "Modules are defined on abundance centred within each person, so they describe",
        "co-movement within people rather than differences between them."
      ),
      paste(
        "Eigengenes (`moduleEigengenes()`) are scored on uncentred abundance in all 45",
        "biopsies; kME (`signedKME()`) is each protein's correlation with its module's",
        "eigengene. `grey` holds proteins in no module."
      ),
      "`icc` is how much of an eigengene's variance lies between people (`lme4::lmer`)."
    ),
    care = c(
      paste(
        "Module colours are names, not ranks, and change when 01_Processing is rebuilt; compare",
        "modules by membership."
      )
    )
  ),
  "01_Processing/07_Module_Characterisation" = list(
    question = paste(
      "What is each module: which gene sets it holds, which proteins lead it, and do its",
      "members interact?"
    ),
    method = c(
      paste(
        "Over-representation of every gene set in each module",
        "(`clusterProfiler::compareCluster(fun = \"enricher\")`), against the measured genes as",
        "universe, BH within module."
      ),
      "Top members: each module's ten highest-kME proteins.",
      paste(
        "STRING v12 edges at combined score 700 or more (`STRINGdb`, files in",
        "`00_Input/STRING`): observed against expected edges per module, with the mapped",
        "measured proteins as background."
      )
    ),
    care = c(
      paste(
        "`top_set` is the lowest p across all collections and a starting point for a name; read",
        "the module's rows in `ora` before naming it."
      ),
      "Module members co-vary by construction, so ORA p-values run small."
    )
  ),
  "02_HRvLR_groups/01_Design_Contrasts" = list(
    question = "Which comparisons does every later model make, and how is that model built?",
    method = c(
      paste(
        "`model.matrix(~ 0 + group + blood)`: six cell means (group by time point) and the",
        "blood index."
      ),
      paste(
        "Nine contrasts (`limma::makeContrasts()`): each group's training (T2 - T1) and acute",
        "(T3 - T2) change, the groups compared at T1, T2 and T3, and the two interactions (one",
        "group's change minus the other's)."
      ),
      paste(
        "Roles: `Training_Interaction` primary, `Acute_Interaction` secondary, the baseline",
        "contrast the floor (the groups come from an outcome, so a T1 difference shows how far",
        "the split separates people before training)."
      ),
      paste(
        "`design.rds` also names the two groups (`arms`), who is in each (`group_of`) and the",
        "traits that defined them; every later step of the grouping reads it."
      )
    ),
    care = character()
  ),
  "02_HRvLR_groups/02_Protein_Contrasts" = list(
    question = "Which proteins change between time points, and which differ between the groups?",
    method = c(
      paste(
        "`fit_limpa()` (`R/helpers.R`): `limpa::dpcDE()` with person as a random block and",
        "limpa's standard errors as weights, then `contrasts.fit()` and `eBayes(robust =",
        "TRUE)`."
      ),
      paste(
        "Each logFC with its 95% interval (`topTable(confint = TRUE)`), BH within contrast,",
        "Storey's q."
      ),
      paste(
        "Checks: the same fit without the blood term (`unadjusted`), null calibration, the",
        "effect detectable at 80% power."
      ),
      paste(
        "`fry` tests one group's p < 0.05 proteins in the other group's contrast",
        "(`fry_concordance`)."
      ),
      "PCA and PERMANOVA (`vegan::adonis2`, 9,999 permutations, time permuted within person)."
    ),
    care = c(
      "Seven or eight people per group: read logFC with its interval, and nominal p as a lead.",
      paste(
        "A protein that changes in one group and not the other is not, on that basis, a group",
        "difference; the interactions test that."
      )
    )
  ),
  "02_HRvLR_groups/03_Set_Contrasts" = list(
    question = "Which gene sets shift in each of the nine contrasts?",
    method = c(
      paste(
        "fgsea (`fgsea::fgsea()`, seeded) ranks each set against the other proteins on the",
        "protein moderated t: does the set sit high or low in the ranking?"
      ),
      paste(
        "fry (`limma::fry()`) tests whether the set's proteins move consistently across people,",
        "on limpa's complete matrix with person blocked."
      ),
      paste(
        "Each method's BH runs within contrast; `fgsea::collapsePathways()` marks the",
        "non-redundant fgsea calls (`main`)."
      )
    ),
    care = c(
      paste(
        "Read the baseline (floor) count first: it is the number of calls the split alone",
        "produces."
      ),
      paste(
        "fgsea treats proteins as exchangeable, which they are not, so its calls run far ahead",
        "of fry's. An fgsea call without a fry call is a ranking result, not evidence that",
        "people moved consistently."
      )
    )
  ),
  "02_HRvLR_groups/04_Module_Contrasts" = list(
    question = "Which co-expression modules shift in each of the nine contrasts?",
    method = c(
      paste(
        "Eigengenes: `lmFit()` with the design from `design.rds`, person blocked",
        "(`duplicateCorrelation()` on the eigengenes), `eBayes(robust = TRUE)`, logFC with its",
        "95% interval."
      ),
      "fry: each module as a protein set on limpa's complete matrix.",
      paste(
        "`membership_contrast`: within each module, Spearman between a member's kME and its",
        "protein t, with a Fisher-z interval."
      )
    ),
    care = c(
      paste(
        "Few modules means BH has little to adjust over; a nominal module result is still one",
        "of several tests per contrast."
      )
    )
  ),
  "02_HRvLR_groups/05_Module_Preservation" = list(
    question = "Does each group's co-expression structure reappear in the other group?",
    method = c(
      paste(
        "Modules built inside each group with the 01_Processing settings",
        "(`WGCNA::blockwiseModules()`)."
      ),
      paste(
        "`WGCNA::modulePreservation()` in both directions, 200 permutations; Zsummary below 2",
        "is none, 2 to 10 weak to moderate, above 10 strong."
      ),
      "Each group module is named by the full-cohort module it overlaps most."
    ),
    care = c(
      paste(
        "Preservation compares correlation structure, not abundance: a module can be preserved",
        "while its level differs between the groups."
      )
    )
  ),
  "02_HRvLR_groups/06_Classification" = list(
    question = "Can a single {noun} tell one group from the other?",
    method = c(
      "Five tasks: the groups compared at T1, T2 and T3, and on the training and acute changes.",
      paste(
        "AUC is the Mann-Whitney W / (n1 n2), above 0.5 when the {unit} is higher in the first",
        "group; Wilcoxon p, BH and Storey q within task{collection}."
      ),
      paste(
        "DeLong 95% interval for every AUC; stratified bootstrap (2,000 resamples) for the",
        "drawn {unit}s (`classify_step()` in `R/classify.R`)."
      )
    ),
    care = c(
      paste(
        "The best AUC among many {noun}s is high even when nothing separates the groups; read",
        "its interval."
      ),
      "With seven or eight people per group the rank-sum p has a floor (`min_p` in `counts`)."
    )
  ),
  "02_HRvLR_groups/07_Trait_Associations" = list(
    question = "Does a protein's, set's or module's link to a trait differ between the groups?",
    method = c(
      "Within each group: Spearman rho with a Fisher-z 95% interval (`spearman_by_row()`).",
      paste(
        "Between groups: limma on one value per person, `feature ~ group * trait`; the",
        "interaction is one group's slope minus the other's, with its 95% interval",
        "(`subject_fit()`)."
      ),
      paste(
        "The traits that defined the groups stay out of the group-difference test. Each group",
        "needs five people with the trait measured."
      )
    ),
    care = c(
      paste(
        "Opposite within-group rho with no group-difference result describes small samples, not",
        "a finding."
      )
    )
  ),
  "03_Cluster_groups/Muscle_size/01_Clustering" = list(
    question = paste(
      "Do people fall into groups by how much their muscle grew, and if two groups are imposed,",
      "who is in each?"
    ),
    method = c(
      paste(
        "PAM on the z-scored % changes over training in fibre area (mixed, type I, type II) and",
        "muscle area, complete cases."
      ),
      paste(
        "K from the gap statistic with the one-SE rule (`cluster::clusGap`, K = 1 to 3, 500",
        "references); stability as bootstrap Jaccard per cluster (`fpc::clusterboot`)."
      ),
      paste(
        "Agreement with HR/LR as the adjusted Rand index (`mclust::adjustedRandIndex`), a",
        "description."
      ),
      paste(
        "When K = 2 is used, `groups.rds` names the clusters Hi and Lo by mean fibre-area",
        "change; `02_Design_Contrasts` reads it."
      )
    ),
    care = c(
      paste(
        "The imposed K = 2 is a decision, not a finding; Hi and Lo differ in size change by",
        "construction."
      )
    )
  ),
  "03_Cluster_groups/Muscle_strength/01_Clustering" = list(
    question = paste(
      "Do people fall into groups by how much their strength rose, and if two groups are",
      "imposed, who is in each?"
    ),
    method = c(
      paste(
        "PAM on the z-scored % changes over training in leg press and leg extension 1RM,",
        "complete cases."
      ),
      paste(
        "K from the gap statistic with the one-SE rule (`cluster::clusGap`, K = 1 to 3, 500",
        "references); stability as bootstrap Jaccard per cluster (`fpc::clusterboot`)."
      ),
      paste(
        "Agreement with HR/LR as the adjusted Rand index (`mclust::adjustedRandIndex`), a",
        "description."
      ),
      paste(
        "When K = 2 is used, `groups.rds` names the clusters Hi and Lo by mean leg-press",
        "change; `02_Design_Contrasts` reads it."
      )
    ),
    care = c(
      paste(
        "The imposed K = 2 is a decision, not a finding; Hi and Lo differ in strength change by",
        "construction."
      )
    )
  ),
  "03_Cluster_groups/Proteome/01_Clustering" = list(
    question = paste(
      "Does the proteome alone, with no phenotype or group label, split these people into",
      "groups?"
    ),
    method = c(
      paste(
        "PAM on the 15 module eigengenes per person, z-scored, at five inputs: T1, T2 and T3",
        "levels, training change (T2 - T1) and acute change (T3 - T2)."
      ),
      paste(
        "K from the gap statistic with the one-SE rule; stability as bootstrap Jaccard;",
        "agreement with HR/LR as the adjusted Rand index."
      ),
      paste(
        "`groups/<input>.rds` holds each input's clusters, C1 and C2 (PAM's numbering, no",
        "direction)."
      )
    ),
    care = c(
      paste(
        "A cluster built from one time point differs in that proteome by construction; later",
        "steps tag every comparison circular, caution or clean (`R/inputs.R`)."
      ),
      "C1 at one input is unrelated to C1 at another."
    )
  ),
  "03_Cluster_groups/Proteome/08_Cluster_Phenotype" = list(
    question = paste(
      "Do people the proteome put in different clusters differ in muscle size or strength?"
    ),
    method = c(
      paste(
        "Per input, kind and trait: `wilcox.test(value ~ cluster)` on the training changes and",
        "the T1 and T2 levels; BH within input."
      ),
      "`hr_lr`: HR and LR counts in each cluster with Fisher's exact p, a description."
    ),
    care = c(
      "At most sixteen people split in two: a large p is not evidence of no difference."
    )
  ),
  "04_Continuous_no_group/01_Trait_Associations" = list(
    question = "Which {noun}s track each size and strength trait, with every person pooled?",
    method = c(
      paste(
        "Sample model: limma on every T1 and T2 biopsy with person blocked; `within` (a",
        "biopsy's trait against that person's mean) and `between` (the person's mean) slopes",
        "with 95% intervals."
      ),
      paste(
        "Change score: one value per person and window; limma slope with its interval, Spearman",
        "rho with a Fisher-z interval, and bicor."
      ),
      paste(
        "Two predictor (`~ growth + strength`) and interaction (`~ growth * strength`) per",
        "change window and 1RM. All in `associate_step()` (`R/associate.R`); BH within each",
        "screen{collection}."
      )
    ),
    care = c(
      "Sixteen people: read every slope with its interval.",
      paste(
        "`within` is the stronger evidence for a link to change, because each person is their",
        "own control."
      )
    )
  ),
  "04_Continuous_no_group/02_Trait_Pathways" = list(
    question = "Which gene sets lead each trait's protein ranking, with every person pooled?",
    method = c(
      paste(
        "fgsea on the limma t from `01_Trait_Associations/01_Proteins` (`change_score`), one",
        "ranking per window and trait, seeded; BH within window and trait."
      )
    ),
    care = c(
      paste(
        "fgsea counts run high (proteins are not exchangeable) and sets overlap, so several",
        "hits can be one signal. No self-contained test (fry) runs here."
      )
    )
  ),
  "04_Continuous_no_group/03_Profile_Clusters" = list(
    question = "Which proteins share a pattern of trait links, with every person pooled?",
    method = c(
      paste(
        "Each protein's limma t across the 42 window-and-trait screens; average-linkage tree on",
        "1 - Pearson r, cut by `dynamicTreeCut::cutreeDynamic(minClusterSize = 20)`."
      ),
      paste(
        "Within-trait groups: the same on each trait's p < 0.05 proteins in `training` and",
        "`baseline_level` (minimum 8)."
      ),
      paste(
        "Each group's most common module and its best over-represented gene set",
        "(`clusterProfiler::enricher()`)."
      )
    ),
    care = c(
      "The 42 screens are not independent; a group can follow a block of related screens.",
      "`top_set` is the best of many sets tried for that group."
    )
  )
)
