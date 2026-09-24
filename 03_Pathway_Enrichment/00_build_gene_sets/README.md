# 00_build_gene_sets

Freezes the gene sets, maps proteins to genes and applies the size filter.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds`, `c_data/cache/goslim_generic.obo` and `c_data/cache/msigdb_2026.1.Hs_Hallmark-Reactome-KEGG_Legacy-GOBP.rds`, each with its `.md5` |
| Writes | `c_data/gene_sets.rds`, `c_data/00_build_gene_sets.xlsx`; the MSigDB cache on the first run |
| Run | `Rscript 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.R` |
| Cost | about 8 s |

`gene_sets.rds` holds `sets` (the 1,378 tested sets as gene symbols), `set_catalog` (every set with
its sizes), `protein_map` (every protein, its symbol and the representative flag) and
`gene_universe` (the 1,900 measured symbols).
