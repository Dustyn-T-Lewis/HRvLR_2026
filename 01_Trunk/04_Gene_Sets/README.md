# 01_Gene_Sets

Freezes the gene sets, maps proteins to genes and applies the size filter.

| | |
|---|---|
| Reads | `01_Trunk/02_Normalisation/c_data/DAList_normalized.rds`, `c_data/cache/goslim_generic.obo` and `c_data/cache/msigdb_2026.1.Hs_Hallmark-Reactome-KEGG_Legacy-GOBP.rds`, each with its `.md5` |
| Writes | `c_data/gene_sets.rds`, `c_data/01_gene_sets.xlsx`; the MSigDB cache on the first run |
| Run | `Rscript 01_Trunk/04_Gene_Sets/a_script/04_gene_sets.R` |

## Method

Four MSigDB collections (Hallmark, KEGG_Legacy, Reactome, GOBP; release 2026.1.Hs) come from a
frozen RDS in `c_data/cache/`, checked against its `.md5`; the first run fetches it with msigdbr. GO
Slim sets come from `goslim_generic.obo` (GO release 2026-07-26): each holds every measured gene
annotated to the slim term or any GO:BP term below it. A set is tested with 15 to 500 source genes
and at least 15 measured. To refresh the OBO, download it into `c_data/cache/` from
<https://current.geneontology.org/ontology/subsets/goslim_generic.obo>, run `md5 -q
goslim_generic.obo > goslim_generic.obo.md5`, and record the release here.

## Figures

None.

## Workbook

| Sheet | Holds |
|---|---|
| `collection_summary` | Sets per collection in the source, sets that qualify, median measured size. |
| `set_catalog` | Every set with its source and measured size and whether it qualifies (15 to 500). |
| `protein_gene_map` | Every protein, its gene symbol, and the representative chosen per symbol. |
| `mapping_summary` | Proteins per mapping status. |
