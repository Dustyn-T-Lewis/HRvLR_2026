# 04 · Network

The module level. Finds proteins that move together, runs the protein level's contrasts, classify
and associate steps on the module eigengenes, then names the modules and asks whether the two arms
share them.

| Step | Runs | Writes |
|---|---|---|
| [`01_Build_Modules`](01_Build_Modules/README.md) | WGCNA on subject-centred abundance | `modules.rds` |
| [`02_Contrasts`](02_Contrasts/README.md) | nine contrasts on eigengenes and with fry, membership against significance | workbook, 5 pages |
| [`03_Classify`](03_Classify/README.md) | AUC and Wilcoxon p on eight tasks, an ROC curve per nominal module | workbook and PDF |
| [`04_Associate`](04_Associate/README.md) | sample-level limma model per trait, change-score Spearman | workbook and PDF |
| [`05_Characterise`](05_Characterise/README.md) | ORA against the stage 03 sets, hubs, STRING enrichment, hub networks | workbook, 5 pages |
| [`06_Preserve`](06_Preserve/README.md) | modules built in each arm, preservation in the other | workbook, 1 page |

## Twelve modules; 285 of 1,900 proteins unassigned

Modules are defined on abundance centred within subject and scored on raw abundance. On raw
abundance subject identity drives the leading components, so modules built there would encode who
a biopsy came from. HR_S28 has one biopsy after outlier removal and centres to zeros, so 44 samples
define the modules and all 45 are scored.

Signed network, biweight midcorrelation (`maxPOutliers = 0.05`), `deepSplit = 2`,
`minModuleSize = 30`, `mergeCutHeight = 0.15`. The soft power is pickSoftThreshold's estimate at
R2 0.85: power 8, R2 0.866, mean connectivity 28.

| Module | Proteins | Subject ICC |
|---|---:|---:|
| turquoise | 363 | 0.00 |
| blue | 182 | 0.17 |
| brown | 170 | 0.00 |
| yellow | 144 | 0.19 |
| green | 134 | 0.00 |
| red | 133 | 0.00 |
| black | 102 | 0.00 |
| pink | 94 | 0.00 |
| magenta | 93 | 0.06 |
| purple | 81 | 0.37 |
| greenyellow | 63 | 0.15 |
| tan | 56 | 0.35 |

## Eight of twelve modules have a set at FDR < 0.05

Enrichment is `compareCluster(fun = "enricher")` against the 1,378 sets of `00_Gene_Sets`,
with the 1,900 measured genes as universe and BH within each module. Hubs are the ten members with
the highest kME. The STRING check is `get_ppi_enrichment()` at combined score ≥ 700 with the
measured proteins STRING maps as background; the expected edge count comes from each member's
degree.

| Module | Top set | FDR | STRING ratio |
|---|---|---:|---:|
| purple | Reactome striated muscle contraction | 7e-27 | 10.1 |
| tan | Reactome respiratory electron transport | 9e-19 | 4.4 |
| pink | Hallmark epithelial-mesenchymal transition | 0.006 | 3.1 |
| red | GO Slim mRNA metabolic process | 1e-4 | 2.8 |
| greenyellow | Reactome rRNA processing | 3e-24 | 2.6 |
| blue | Hallmark oxidative phosphorylation | 1e-15 | 1.9 |
| brown | Hallmark glycolysis | 0.006 | 1.8 |
| magenta | Reactome eukaryotic translation initiation | 3e-4 | 1.7 |

Turquoise, yellow, green and black have no set at FDR < 0.05. All twelve modules share more STRING
edges than their degrees predict (1.1 to 10.1 times; FDR < 0.05). Among the 25 highest-kME members,
greenyellow has 209 edges, tan 202 and purple 167; turquoise, green and black have two to four.

## 35 of 38 arm modules reach at least moderate preservation

Each arm gets its own network with the settings above. HR gives 21 modules at power 16 (the WGCNA
FAQ fallback for 20 samples, since no power cleared 0.85) and LR 17 at power 9.
`modulePreservation()` tests each arm's modules in the other, 200 permutations. Zsummary above 10
is strong preservation, 2 to 10 moderate, below 2 none (Langfelder et al. 2011).

| Direction | Strong | Moderate | None |
|---|---:|---:|---:|
| HR modules in LR | 3 | 17 | 1 |
| LR modules in HR | 4 | 11 | 2 |

Turquoise is strongly preserved both ways (Zsummary 28 and 31), as is the muscle contraction
module (HR purple 15.3 in LR; LR greenyellow, which maps to full-cohort purple, 16.2 in HR). Three
modules of 39 to 54 proteins fall below 2: HR lightyellow, and LR lightcyan and tan. Each arm
network rests on 20 to 24 samples.

## No eigengene clears BH; fry finds red in Acute_LR

The eigengene test is `lmFit()` with the protein design and subject block, the correlation
re-estimated on the eigengenes (0.095), `eBayes(robust = TRUE)` and BH within contrast. fry tests
each module as a protein set on the imputed matrix with a correlation estimated on that matrix
(0.176); its FDR runs over the twelve modules within each contrast.

No eigengene clears BH in any contrast, task, trait or outcome. fry finds one module: red rises
in Acute_LR (FDR 0.045).

| Test | Tests | Expected at p < 0.05 | Nominal | BH < 0.05 |
|---|---:|---:|---:|---:|
| eigengene contrasts | 108 | 5.4 | 5 | 0 |
| fry contrasts | 108 | 5.4 | 7 | 1 |
| classification | 96 | 4.8 | 6 | 0 |
| sample model | 216 | 10.8 | 9 | 0 |
| change score | 720 | 36 | 25 | 0 |

In the sample model purple, the striated-muscle module, rises within-person with type I fibre
share (p = 0.034), the direction the protein-level positive control shows.

Membership against significance is the Spearman correlation, within each module, between a
member's kME and its protein-level moderated t. On the primary contrast blue (rho 0.32) and green
(0.32) run positive and pink negative (−0.32); on the acute interaction green runs negative
(−0.51). Members share a module, so these correlations are descriptive.
