# 02_characterise_modules

Names each module by its enriched sets, hubs and STRING edges.

| | |
|---|---|
| Reads | `01_build_modules/c_data/modules.rds`, `03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds`, `00_Input/downloads/9606.protein.links.v12.0.min700.txt.gz`, `9606.protein.aliases.v12.0.txt.gz`, `9606.protein.info.v12.0.txt.gz` |
| Writes | `c_data/02_characterise_modules.xlsx`, `b_reports/02_characterise_modules_figures.pdf` |
| Run | `Rscript 04_Network/02_characterise_modules/a_script/02_characterise_modules.R` |
| Cost | about 15 s with the STRING files on disk |

The figure PDF holds the enrichment dot plot, observed against expected STRING edges, and three
pages of hub networks: each module's 25 highest-kME members and the STRING edges among them, four
modules to a page. The seeded layout carries no meaning. The `labels` sheet holds each module's top
set, top GO Slim term and hubs.
