# 02_characterise_modules

What each module is.

| | |
|---|---|
| **Reads** | `modules.rds`, `gene_sets.rds`, STRING v12 files in `00_Input/downloads/` |
| **Writes** | `module_characterisation.rds`, `02_characterise_modules.xlsx`, 2 figures |

**Enrichment.** `clusterProfiler::compareCluster(fun = "enricher")` against the 1,378 sets of
`03_Pathway_Enrichment/00_build_gene_sets`, with the 1,900 measured genes as universe. BH within
each module.

**Hubs.** The ten members with the highest kME.

**STRING.** `STRINGdb$get_ppi_enrichment()` at combined score ≥ 700, background set to the
measured proteins STRING maps. Its expected edge count comes from each member's degree, so a
module of well-studied proteins is not rewarded for being well studied. Nearly any co-expression
module passes this test; the observed-over-expected ratio carries the information.

| Module | Top set | STRING ratio |
|---|---|---:|
| pink | Reactome striated muscle contraction | 10.1 |
| tan | Reactome respiratory electron transport | 4.2 |
| magenta | Reactome GCN2 response to amino acid deficiency | 2.6 |
| brown | GO Slim cytoskeleton organization | 2.4 |
| yellow | Reactome aerobic respiration and electron transport | 2.0 |
| greenyellow | Hallmark epithelial-mesenchymal transition | 2.0 |

The full table, with every module's hubs, is the workbook's `labels` sheet.
