# 07_Module_Characterisation

What is each module about?

## The idea

`06_Module_Building` returns modules as colours and member lists. A colour says nothing about
biology, so each module is described three ways, each answering a different question:

| Evidence | Question | Method |
|---|---|---|
| Enriched gene sets | Do the members share known functions more than chance? | Over-representation analysis (ORA) |
| Hubs | Which proteins track the module most closely? | Highest kME |
| STRING edges | Do the members physically or functionally interact more than chance? | STRING PPI enrichment |

ORA asks, for each gene set: of the module's members, how many belong to the set, and how
surprising is that count given the set's size and the module's size? The hypergeometric test gives
the p-value. The background matters. Every protein here was measured in muscle, so testing against
the whole genome would make any module look enriched for muscle. The background is the gene
universe from `04_Gene_Sets`: the measured genes, one per symbol. BH adjusts p within each module.

STRING is a database of known and predicted protein-protein associations, each with a combined
confidence score from 0 to 1000. At a score of 700 or more (high confidence), the test counts the
edges among a module's members and compares the count with the number expected from how connected
each member is in the background. A module with far more edges than expected is held together by
known biology as well as by correlation. The background is the measured proteins STRING can map,
not the whole proteome.

Hubs are the module's highest-kME proteins: those whose abundance follows the eigengene most
closely within a person. They are the module's most typical members, not necessarily its most
important ones.

## What the code does

| Step | Function | Detail |
|---|---|---|
| ORA | `clusterProfiler::compareCluster(fun = "enricher")` | Each module's member genes (grey excluded) against the qualifying sets, universe `gene_universe`, `minGSSize = 15`, `maxGSSize = 500`, no p or q cut-off, so every tested set is kept. BH within module. |
| Hubs | `dplyr::slice_max(kme, n = 10)` | The ten highest-kME members of each module. |
| STRING map | `STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700)`, `$map()` | Reads the local files in `00_Input/downloads`. Every module protein, grey included, is mapped by UniProt accession; unmapped proteins drop out. `$set_background()` sets the mapped proteins as background. |
| STRING enrichment | `$get_ppi_enrichment()` | Per module: observed edges, expected edges (`lambda`), their ratio and a p; BH across modules. |
| Labels | `dplyr::slice_min(p)` | Per module: the set with the lowest ORA p, the GO Slim term with the lowest p, the hubs and the STRING ratio. |
| Hub networks | `$get_interactions()`, `tidygraph`, `ggraph` | Each module's 25 highest-kME members and the STRING edges among them at score 700 or more. Force-directed layout, seeded. |

## How to read the outputs

- `labels` is a starting point for a name, not a name. `top_set` is the lowest p across all
  collections, so a large GOBP term can win over a precise Hallmark set; read the module's rows in
  `ora` before naming it.
- The hypergeometric test treats genes as independent draws. Module members were chosen because
  they co-vary, and members of one pathway co-vary too, so ORA p-values within modules run small.
- `string` ratio is observed over expected edges. A small module with few expected edges can show a
  large ratio from a handful of edges; read `n_edges` beside it.
- The hub-network layout carries no meaning: distance and position on the page are arbitrary.
  Node fill is kME, node size the number of edges among the drawn hubs, edge width the STRING
  score.
- Module colours change when the trunk is rebuilt, so this workbook belongs to one build of
  `modules.rds`.

## Files

| | |
|---|---|
| Reads | `01_Trunk/06_Module_Building/c_data/modules.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `00_Input/downloads/9606.protein.links.v12.0.min700.txt.gz`, `00_Input/downloads/9606.protein.aliases.v12.0.txt.gz`, `00_Input/downloads/9606.protein.info.v12.0.txt.gz` |
| Writes | `c_data/module_characterisation.xlsx` (read by `02_HRvLR_groups/Figures/F2_Programmes`), `b_reports/module_characterisation_report.pdf` |
| Run | `Rscript 01_Trunk/07_Module_Characterisation/a_script/module_characterisation.R` |

The STRING files are not in git. `00_Input/README.md` has the download command; the script stops
and names the missing files if they are absent.

```
07_Module_Characterisation/
├── a_script/module_characterisation.R
├── b_reports/module_characterisation_report.pdf
└── c_data/module_characterisation.xlsx
```

## Workbook: `module_characterisation.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close it.

| Sheet | One row is | Key columns |
|---|---|---|
| `labels` | a module | `n_proteins`, `top_set`, `top_set_fdr`, `top_go_slim`, `go_slim_fdr`, `hubs`, `string_ratio`, `string_fdr` |
| `ora` | a module and gene set | `set_id`, `collection`, `pathway`, `gene_ratio`, `bg_ratio`, `n_overlap`, `p`, `fdr` |
| `hubs` | a hub protein | `module`, `uniprot_id`, `gene`, `kme` |
| `string` | a module | `n_mapped`, `n_edges`, `n_expected`, `ratio`, `p`, `fdr` |
| `hub_edges` | a STRING edge between two drawn hubs | `module`, `from_gene`, `to_gene`, `combined_score` |

## Report: `module_characterisation_report.pdf`

1. Enrichment dot plot: the top three sets per module among those at BH < 0.05; dot size is gene
   ratio, colour BH. Modules with no set at BH < 0.05 do not appear.
2. STRING edges per module: observed edges against expected, on a log scale, labelled with the
   ratio.
3. Hub networks, four modules to a page, as many pages as the modules need.
