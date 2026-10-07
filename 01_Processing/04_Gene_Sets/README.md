# 04_Gene_Sets

Which groups of proteins count as one biological process, decided once for every arc?

| | |
|---|---|
| Reads | `01_Processing/02_Normalisation/c_data/DAList_normalized.rds`, `00_Input/msigdb_2026.1.Hs.rds`, `00_Input/goslim_generic.obo` |
| Writes | `c_data/gene_sets.rds`, `c_data/gene_sets.xlsx` |
| Run | `Rscript 01_Processing/04_Gene_Sets/a_script/gene_sets.R` |

## Key result

- 1,378 sets qualify for testing: Hallmark 33 of 50, Reactome 328 of 1,839, KEGG_Legacy 62 of 186, GOBP 897 of 7,538, GO_Slim 58 of 71.
- All 1,900 proteins map to one distinct gene symbol each.

## Method

- MSigDB 2026.1.Hs (Hallmark, Reactome, KEGG Legacy, GO:BP) from the snapshot in `00_Input`, plus GO Slim terms built from the GO:BP hierarchy (`GSEABase`, `GO.db`).
- One protein per gene symbol: the one observed in the most samples. Proteins with no symbol or several are left out.
- A set is tested when it holds 15 to 500 genes in the source and at least 15 measured here (GO Slim: 15 to 500 measured).

## Read with care

- Overlapping sets are kept, so one signal can appear in several neighbouring sets.

## Sheets

`gene_sets.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `collection_summary` | Sets per collection in the source, sets that qualify, median measured size. |
| `set_catalog` | Every set with its source and measured size and whether it qualifies (15 to 500). |
| `protein_gene_map` | Every protein, its gene symbol, and the representative chosen per symbol. |
| `mapping_summary` | Proteins per mapping status. |
