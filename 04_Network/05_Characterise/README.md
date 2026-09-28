# 05_Characterise

Names each module by its enriched sets, hubs and STRING edges.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds`, `03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds`, `00_Input/downloads/9606.protein.links.v12.0.min700.txt.gz`, `9606.protein.aliases.v12.0.txt.gz`, `9606.protein.info.v12.0.txt.gz` |
| Writes | `c_data/05_characterise.xlsx`, `b_reports/05_characterise_figures.pdf` |
| Run | `Rscript 04_Network/05_Characterise/a_script/05_characterise.R` |

## Method

`clusterProfiler::compareCluster(fun = "enricher")` against the `00_Gene_Sets` sets with the
measured genes as universe, BH within module. Hubs are the ten highest-kME members. `STRINGdb`
`get_ppi_enrichment()` at combined score 700 or more, with the measured proteins STRING maps as
background. The STRING files must be on disk; `00_Input/README.md` has the download.

## Figures

The enrichment dot plot, observed against expected STRING edges, and three pages of hub networks:
each module's 25 highest-kME members and their STRING edges, four modules to a page. The seeded
layout carries no meaning.

## Workbook

| Sheet | Holds |
|---|---|
| `labels` | One row per module: size, top set, top GO Slim term, hubs, STRING enrichment. |
| `ora` | clusterProfiler ORA of every module against the 00_Gene_Sets sets, BH within module. |
| `hubs` | The ten highest-kME members of each module. |
| `string` | STRINGdb PPI enrichment per module: observed and expected edges at score >= 700. |
| `hub_edges` | STRING edges among each module's 25 highest-kME members, as drawn. |
