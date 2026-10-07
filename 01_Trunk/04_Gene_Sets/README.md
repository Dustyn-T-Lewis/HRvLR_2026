# 04_Gene_Sets

Which groups of proteins count as one biological process, decided once for every arc?

## The idea

One protein moving a little is hard to tell from noise. Thirty proteins that run the same process
all moving a little the same way is a stronger signal. A gene set is such a group: a list of genes
that a curated database assigns to one pathway, process or cell state. Set-level tests pool the
members, so they can see a shift no single protein clears.

The sets come from MSigDB, which collects them from several sources:

| Collection | What a set is |
|---|---|
| Hallmark | a well-defined biological state or process, refined from many overlapping sets |
| Reactome | a pathway in the Reactome database |
| KEGG_Legacy | a pathway in the older, freely licensed KEGG release |
| GOBP | a Gene Ontology biological-process term |
| GO_Slim | a broad GO term from the GO Consortium's generic slim, holding every gene under it in the GO:BP hierarchy |

Three choices are made here so that no later step makes them differently:

1. One release. Set membership changes between MSigDB releases, so release 2026.1.Hs is frozen in
   `c_data/cache/` with an md5 checksum. The GO Slim file is frozen the same way.
2. One protein per gene. Sets list gene symbols; the matrix holds proteins. A protein with no
   symbol or several is left out. Where several proteins share a symbol, the one observed in the
   most samples represents it (ties go to more peptide sequences, then to the accession). This is
   decided before any fold change exists.
3. One size rule. A set is tested when it holds 15 to 500 genes in the source and at least 15
   genes measured here. Below 15 a set is a handful of proteins; above 500 it says little about
   any one process. GO Slim terms are broad by design, so their 15-to-500 rule reads the measured
   size. No set is dropped for overlapping another or for its name.

The measured genes, one per symbol, form the gene universe. `07_Module_Characterisation` uses it
as the background for its enrichment test, so a module is compared with the muscle proteome that
was measured rather than with the whole genome.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Fetch (first run only) | `msigdbr::msigdbr()` | Hallmark (`H`), Reactome (`C2`, `CP:REACTOME`), KEGG_Legacy (`C2`, `CP:KEGG_LEGACY`), GOBP (`C5`, `GO:BP`). Stops if msigdbr returns a release other than 2026.1.Hs. Writes the snapshot and its `.md5`. |
| Verify | `tools::md5sum()` | Every run checks the snapshot and the GO Slim file against their `.md5`; a mismatch stops the run. Later runs need no network. |
| Map | `dplyr::distinct()` | `protein_map`: each protein's symbol, mapping status and the representative per symbol; `label` adds the accession where two proteins share a symbol. |
| Size filter | `lengths()` | Source and measured size per set; `qualifies` marks the sets kept. |
| GO Slim | `GSEABase::getOBOCollection()`, `GO.db::GOBPOFFSPRING` | Each slim term gathers the measured genes annotated to it or any GO:BP term beneath it, from the full frozen GOBP membership, so a gene is not lost when its own GOBP set fails the size filter. |

## How to read the outputs

- `collection_summary` counts sets in the source and sets that qualify, with the median measured
  size of those that qualify.
- In `set_catalog`, `source_size` against `measured_size` shows how much of a set this experiment
  sees. A set that qualifies with a small share measured rests on the proteins that were detected,
  which lean towards the abundant ones.
- `mapping_status` is `representative`, `duplicate_gene`, `missing_symbol` or `multiple_symbols`.
  Only representatives enter the sets and the universe.
- The script reads `DAList_normalized.rds` for annotation and detection counts only, so the sets do
  not depend on limpa.

## Files

| | |
|---|---|
| Reads | `01_Trunk/02_Normalisation/c_data/DAList_normalized.rds`, `c_data/cache/goslim_generic.obo` (with `.md5`), `c_data/cache/msigdb_2026.1.Hs_Hallmark-Reactome-KEGG_Legacy-GOBP.rds` (with `.md5`) |
| Writes | `c_data/gene_sets.rds`, `c_data/gene_sets.xlsx`; the MSigDB snapshot on the first run |
| Run | `Rscript 01_Trunk/04_Gene_Sets/a_script/gene_sets.R` |

`gene_sets.rds` holds `sets` (the qualifying sets as gene vectors, named by `set_id`),
`set_catalog`, `protein_map` and `gene_universe`. `05_Set_Scores`, `07_Module_Characterisation`
and every set-level step in the three arcs read it.

```
04_Gene_Sets/
├── a_script/gene_sets.R
└── c_data/
    ├── gene_sets.xlsx
    ├── gene_sets.rds
    └── cache/
        ├── goslim_generic.obo
        ├── goslim_generic.obo.md5
        ├── msigdb_2026.1.Hs_Hallmark-Reactome-KEGG_Legacy-GOBP.rds
        └── msigdb_2026.1.Hs_Hallmark-Reactome-KEGG_Legacy-GOBP.rds.md5
```

`set_id` is `collection|name`, for example `Hallmark|HALLMARK_MYOGENESIS`. A GO Slim name is the
term name in capitals with each run of other characters replaced by an underscore.

## Restoring or refreshing the GO Slim file

The cached `goslim_generic.obo` is GO release 2026-07-26 (its `data-version` line). To refresh it,
download it into `c_data/cache/` from
<https://current.geneontology.org/ontology/subsets/goslim_generic.obo>, then run
`md5 -q goslim_generic.obo > goslim_generic.obo.md5` and record the new release here. Deleting the
MSigDB snapshot and its `.md5` makes the next run fetch it again.

## Workbook: `gene_sets.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close it.

| Sheet | One row is | Key columns |
|---|---|---|
| `collection_summary` | a collection | `n_in_source`, `n_qualifying`, `median_measured` |
| `set_catalog` | a set in the source | `set_id`, `collection`, `pathway`, `source_id`, `source_size`, `measured_size`, `qualifies` |
| `protein_gene_map` | a protein | `gene`, `mapping_status`, `selected`, `label`, `n_obs`, `n_seq` |
| `mapping_summary` | a mapping status | `n_proteins` |

No report PDF.
