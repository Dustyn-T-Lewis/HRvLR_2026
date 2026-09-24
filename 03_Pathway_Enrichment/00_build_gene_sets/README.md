# 00_build_gene_sets

Freezes the gene sets, maps proteins to genes, and applies the size filter. Every later step reads
this one list.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds`, `c_data/cache/goslim_generic.obo` |
| Writes | `gene_sets.rds`, `00_build_gene_sets.xlsx`, `c_data/cache/` |

## Snapshot

The first run fetches four MSigDB collections (release 2026.1.Hs) and writes an RDS with an md5.
Later runs verify the checksum and need no network. Restore the RDS and its `.md5` together; move
both aside to rebuild. `goslim_generic.obo` (GO release 2026-07-26) sits beside it with its own
`.md5`; the script does not fetch it. To replace it, download it from
<https://current.geneontology.org/ontology/subsets/goslim_generic.obo> into the cache and rewrite
the checksum (`00_Input/PROVENANCE.md`).

## Collections

| Collection | In source | Tested | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 33 | 34 |
| KEGG_Legacy | 186 | 62 | 23.5 |
| Reactome | 1,839 | 328 | 33 |
| GO:BP | 7,538 | 897 | 25 |
| GO Slim | 71 | 58 | 97 |
| Total | 9,684 | 1,378 | |

A set is tested with 15 to 500 source genes and at least 15 measured. GO Slim sets are built here:
each holds every measured gene annotated to the slim term or any GO:BP term below it, taken from
the full membership, and the size rule reads measured size.

## One protein per gene

All 1,900 proteins carry one distinct symbol, so every protein represents its gene. The rule for
a shared symbol (keep the protein seen in the most samples) is in place for a future matrix.
`protein_gene_map` records each decision.

```r
gs <- readRDS("03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds")
gs$sets          # 1,378 tested sets, measured gene symbols
gs$set_catalog   # every set, with sizes and whether it is tested
gs$protein_map   # representative protein per gene, and plot labels
gs$gene_universe # 1,900 measured symbols
```
