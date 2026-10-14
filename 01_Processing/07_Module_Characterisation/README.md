# 07_Module_Characterisation

What is each module: which gene sets it holds, which proteins lead it, and do its members interact?

| | |
|---|---|
| Reads | `01_Processing/06_Module_Building/c_data/modules.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/07_Module_Characterisation/references/STRING/9606.protein.links.v12.0.min700.txt.gz`, `01_Processing/07_Module_Characterisation/references/STRING/9606.protein.aliases.v12.0.txt.gz`, `01_Processing/07_Module_Characterisation/references/STRING/9606.protein.info.v12.0.txt.gz` |
| Writes | `c_data/module_characterisation.xlsx`, `b_reports/module_characterisation_report.pdf` |
| Run | `Rscript 01_Processing/07_Module_Characterisation/a_script/module_characterisation.R` |

## Key result

- 11 of 15 modules have a gene set over-represented at FDR < 0.05; 14 have more STRING edges than expected at FDR < 0.05.
- Strongest: `pink`, REACTOME_STRIATED_MUSCLE_CONTRACTION (FDR 3e-32); `red`, REACTOME_EUKARYOTIC_TRANSLATION_INITIATION (FDR 9e-32); `midnightblue`, REACTOME_RESPIRATORY_ELECTRON_TRANSPORT (FDR 1e-22).

## Method

- Over-representation of every gene set in each module (`clusterProfiler::compareCluster(fun = "enricher")`), against the measured genes as universe, BH within module.
- Top members: each module's ten highest-kME proteins.
- STRING v12 edges at combined score 700 or more (`STRINGdb`, files in `01_Processing/07_Module_Characterisation/references/STRING`): observed against expected edges per module, with the mapped measured proteins as background.

## Read with care

- `top_set` is the lowest p across all collections and a starting point for a name; read the module's rows in `ora` before naming it.
- Module members co-vary by construction, so ORA p-values run small.

## Sheets

`module_characterisation.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `labels` | One row per module: size, top set, top GO Slim term, top members, STRING enrichment. |
| `ora` | clusterProfiler ORA of every module against the 04_Gene_Sets sets, BH within module. |
| `top_members` | The ten highest-kME members of each module. |
| `string` | STRINGdb PPI enrichment per module: observed and expected edges at score >= 700. |
| `network_edges` | STRING edges among each module's 25 highest-kME members, as drawn. |
