# 00_Gene_Sets

Freezes the gene sets, maps proteins to genes and applies the size filter.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds`, `c_data/cache/goslim_generic.obo` and `c_data/cache/msigdb_2026.1.Hs_Hallmark-Reactome-KEGG_Legacy-GOBP.rds`, each with its `.md5` |
| Writes | `c_data/gene_sets.rds`, `c_data/00_gene_sets.xlsx`; the MSigDB cache on the first run |
| Run | `Rscript 03_Pathway_Enrichment/00_Gene_Sets/a_script/00_gene_sets.R` |
| Cost | about 8 s |

`gene_sets.rds` holds `sets` (the 1,378 tested sets as gene symbols), `set_catalog` (every set with
its sizes), `protein_map` (every protein, its symbol and the representative flag) and
`gene_universe` (the 1,900 measured symbols).

## GO Slim generic, GO release 2026-07-26

`goslim_generic.obo` is the Gene Ontology Consortium's species-neutral slim, downloaded 2026-08-19
from <https://current.geneontology.org/ontology/subsets/goslim_generic.obo>. It holds 140 terms (72
biological process, 28 cellular component, 40 molecular function); this step builds 71 sets from
the biological-process terms. The script stops if the file's checksum does not match its `.md5`.

To refresh, download into `c_data/cache/`, then run `md5 -q goslim_generic.obo > goslim_generic.obo.md5`
and record the new release here. Each release can change term membership and every result that
reads it.
