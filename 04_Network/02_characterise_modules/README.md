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
measured proteins STRING maps. The expected edge count comes from each member's degree, so
well-studied proteins do not inflate a module's enrichment. Nearly any co-expression
module passes; compare modules by the observed-over-expected ratio.

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

Turquoise, yellow, green and black have no set at FDR < 0.05. The workbook's `labels` sheet holds
every module's top set, top GO Slim term and hubs.
