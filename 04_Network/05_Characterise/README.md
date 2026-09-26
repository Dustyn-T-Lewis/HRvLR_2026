# 05_Characterise

Names each module by its enriched sets, hubs and STRING edges.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds`, `03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds`, `00_Input/downloads/9606.protein.links.v12.0.min700.txt.gz`, `9606.protein.aliases.v12.0.txt.gz`, `9606.protein.info.v12.0.txt.gz` |
| Writes | `c_data/05_characterise.xlsx`, `b_reports/05_characterise_figures.pdf` |
| Run | `Rscript 04_Network/05_Characterise/a_script/05_characterise.R` |
| Cost | about 15 s with the STRING files on disk |

The figure PDF holds the enrichment dot plot, observed against expected STRING edges, and three
pages of hub networks: each module's 25 highest-kME members and the STRING edges among them, four
modules to a page. The seeded layout carries no meaning. The `labels` sheet holds each module's top
set, top GO Slim term and hubs.
