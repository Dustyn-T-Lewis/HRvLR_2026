# 00_Setup

Study inputs and the small helpers every walkthrough sources. Nothing here runs on its own.

| Folder | Holds |
|---|---|
| [`Input`](Input/README.md) | DIA-NN protein matrix, sample sheet, phenotype |
| [`Shared_functions`](Shared_functions/README.md) | helpers used by more than one walkthrough: workbook writing, limpa fit, Storey q, relabelling a grouping, report heatmaps, figure theme and panels, the Proteome input stacker |
| `Manuscript` | the collaborator's manuscript draft; local only, not tracked |

Reference snapshots (MSigDB, GO Slim, STRING, Human Protein Atlas) live in the `references/` folder of
the processing step that reads them. `renv/`, `renv.lock` and `.Rprofile` sit at the project root;
run every command from there.
