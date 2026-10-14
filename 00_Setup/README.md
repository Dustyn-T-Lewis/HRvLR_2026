# Setup and shared code

This folder holds study inputs and code reused by multiple stages. It is not an executable stage: open the numbered processing or analysis substage to run an analysis.

| Folder | Purpose |
|------------------------------------|------------------------------------|
| [Input](Input/README.md) | DIA-NN protein quantities, sample metadata and participant outcomes |
| [Shared_functions](Shared_functions/README.md) | Statistical, reporting and figure helpers |
| [Shared_analyses](Shared_analyses/README.md) | Complete analyses reused across groupings or feature levels |

A substage's `a_script/` file is its entry point. It loads packages, explicitly sources shared code if needed, and supplies the stage-specific inputs and output location. A `run_*()` call means that the analysis implementation lives in `Shared_analyses/`, not in a small helper. These stages are explicitly shared, rather than fully self-contained.

Stage-specific reference snapshots remain in the owning processing stage's `references/`. The Proteome grouping's input helper (`validity()`, `stack_inputs()`) is `Shared_functions/proteome_inputs.R`.

`renv/`, `renv.lock` and `.Rprofile` remain at the repository root to use renv's standard startup and restore layout. Open the root project and run commands from the repository root.