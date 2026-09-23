# 03 · Pathway Enrichment

Tests whether proteins that work together moved together, then scores every sample on every set.

| Step | Runs | Writes |
|---|---|---|
| [`00_build_gene_sets`](00_build_gene_sets/README.md) | frozen MSigDB snapshot, protein-to-gene map, size filter, GO Slim sets | `gene_sets.rds` |
| [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) | fgsea and fry per contrast, `collapsePathways` | `set_tests.rds`, dot plots, hit matrices |
| [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) | protein volcanoes with pathway rings | 8 volcanoes |
| [`03_enrich_scatter_fgsea`](03_enrich_scatter_fgsea/README.md) | HR NES against LR NES, training and acute | 4 composites |
| [`04_run_singscore`](04_run_singscore/README.md) | per-sample set scores | `singscore.rds`, 2 figures |
| [`05_classify_and_associate_sets`](05_classify_and_associate_sets/README.md) | set classification and phenotype association | `set_results.rds`, ROC, association, chance and hit-matrix figures |

```sh
for s in 00_build_gene_sets 01_run_fgsea_and_fry 02_enrich_volcano_fgsea \
         03_enrich_scatter_fgsea 04_run_singscore 05_classify_and_associate_sets; do
  Rscript 03_Pathway_Enrichment/$s/a_script/$s.R
done
```

About two minutes in total. Each substage writes every figure as PNG and PDF and bundles them into
one `<substage>_figures.pdf`. Titles name the figure, subtitles give method and counts, captions
state the encodings and the source table. Findings live here, not on the figures.

## Methods

**fgsea** is competitive: does a set sit at one end of the protein ranking? **fry** is
self-contained: did the set move at all under the fitted design? Both run on the same 1,378 sets.
fgsea assumes proteins are exchangeable, and co-regulated sets break that; fry does not assume it,
and takes the subject block and the within-subject correlation. fry reads the imputed matrix,
because it cannot take a missing value.

**singscore** gives one rank-based score per set per sample, on the imputed matrix. It carries no
p-value, never sees a contrast, and a sample's score does not change with the cohort.

**Order: pool, test, then collapse.** All sets are tested and each method's own BH corrects within
each contrast. `collapsePathways` then marks non-redundant fgsea hits in the `main` column and
deletes nothing. No set is excluded by name.

## Results

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| Training_Interaction *(primary)* | 27 | 13 | 0 |
| Acute_Interaction *(secondary)* | 110 | 44 | 0 |
| Training_HR | 68 | 23 | 0 |
| Training_LR | 121 | 18 | 0 |
| Acute_HR | 279 | 83 | 0 |
| Acute_LR | 42 | 18 | 0 |
| Trained_HRvLR | 0 | 0 | 0 |
| Acute_HRvLR | 175 | 58 | 0 |
| Baseline_HRvLR *(floor)* | 60 | 26 | 0 |

fry finds nothing in any contrast. fgsea calls 60 sets on the floor, more than on the primary
contrast, so its counts here measure its gene-permutation null rather than biology. None of the
fgsea lists is a finding.

**Classification by set score**, nominal hits over chance, per collection:

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| training, HR | 0.59 | 1.29 | 0.79 | 0.60 | 1.72 |
| training, LR | 0.00 | 0.65 | 0.12 | 0.56 | 0.34 |
| acute, HR | 2.35 | 1.94 | 1.34 | 1.71 | 1.38 |
| acute, LR | 2.94 | 4.19 | 1.77 | 2.23 | 2.07 |
| HR vs LR at T1 *(floor)* | 0.59 | 0.00 | 0.12 | 0.24 | 0.00 |
| HR vs LR at T2 | 0.00 | 0.00 | 0.06 | 0.29 | 0.00 |
| HR vs LR, training change | 1.18 | 1.94 | 0.30 | 0.47 | 1.03 |
| HR vs LR, acute change | 1.18 | 0.97 | 0.37 | 0.98 | 0.34 |

The acute bout clears chance in every collection in both arms; training does not; no between-arm
task does. No set survives BH in any task. With 7 or 8 pairs, the paired Wilcoxon cannot go below
p = 0.0078, so BH over hundreds of sets is out of reach by construction.

**Association.** Four set-outcome pairs survive BH within collection: peroxisomal protein import
(Reactome) and peroxisome (KEGG) against `volume_load` over training, GO Slim DNA replication
against baseline `d_1rm_ext`, and Reactome basal-body anchoring against `comp_hypertrophy` over the
acute bout. Baseline level against `d_1rm_ext` runs at 3.4 times chance across collections.

**Concordance.** HR and LR NES correlate at rho 0.34 over all sets for training and 0.35 for the
acute bout. The overlap is weaker than in BFR's two training arms (0.83).
