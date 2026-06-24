# Repository cleanup and alignment

Updated: 2026-06-23

## Goal

Bring `A_HRvLR_2026` to the same clean state as the `A_CvH_2026` and
`A_YvO_2026` reference pipelines, so a new collaborator can navigate it without
a guide. This is consolidation, not a restructure: the repo already uses the
shared convention. The work removes stale duplicates, gathers reused helpers
into `functions/` subdirectories, tidies outputs, and scrubs machine-written
comment noise.

## Conventions confirmed

- Numbered, self-contained stages: `00_input` to `05_Clustering`.
- Each analysis stage and figure directory keeps the triad `a_script/` (code),
  `b_reports/` (rendered PDF/PNG), `c_data/` (csv/xlsx/rds outputs).
- `a_script/` stays singular, matching the reference projects and every existing
  directory.
- Reused code lives in a `functions/` subdirectory that scripts source. Helpers
  used by only one script stay inline; pulling them out would add indirection
  without a second caller.

## Current state

The pipeline runs `01_Filtering -> 02_Normalization -> 03_DEP -> 04_Figures`,
with `05_Clustering` as the newest stage. Stale artifacts have accumulated
alongside the live ones:

| Stale item | Superseded by |
| --- | --- |
| `01_normalization/` (only `.DS_Store`) | `01_Filtering/`, `02_Normalization/` |
| `02_Imputation/` (only `.DS_Store`) | `02_Normalization/imputation/` |
| `03_DEP/{a_script,b_reports,c_data}` (older copy) | `03_DEP/a_non_imputed/` |
| `docs/archive/` (retired F04 to F06) | live `04_Figures/F01` to `F03` |

Reused code today sits in `04_Figures/shared/` (`style.R`, `pathway_utils.R`,
sourced by about two dozen panels) and `05_Clustering/_shared_inputs.R` (sourced
by four scripts). Several F02 panels carry many helper definitions inline; for
example `panel_I.R` defines the venn builders that `panel_J.R` reuses.

## Plan

### Phase 1: remove dead weight

Each target was confirmed unreferenced by any active script.

- Delete `01_normalization/` and `02_Imputation/`.
- Delete the root `03_DEP/{a_script,b_reports,c_data}` copy, keeping
  `a_non_imputed/` and `b_imputed/`.
- Delete `docs/archive/`.
- Remove stray `.DS_Store` files and confirm `.gitignore` covers them.
- Flag `00_input/HPA_skeletal_muscle_annotations.tsv` (3.4 MB, unread) for the
  maintainer to confirm before deletion; filtering reads `HPA_annotations.tsv`.

### Phase 2: gather helpers into `functions/`

- Rename `04_Figures/shared/` to `04_Figures/functions/` and update the
  `source()` paths in every figure that reads it.
- Move `05_Clustering/_shared_inputs.R` to `05_Clustering/functions/inputs.R`
  and update the four callers.
- Pull reused inline helpers out of cluttered panels into a per-figure
  `a_script/functions/` file, so the panel scripts read as a sequence of plot
  steps. The venn builders shared by F02 `panel_I` and `panel_J` are the first
  case.
- Leave the single-use helpers in stages 01 to 03 inline.

### Phase 3: tidy outputs and names

- Keep `b_reports/` to current renders; supplementary panels go under `supp/`,
  the pattern F02 already adopted.
- Keep each stage's outputs and its supplementary `.xlsx` in `c_data/`; drop
  superseded `audit_panel_*.csv` only where a current file replaces them.
- Use one panel and output naming pattern across F01, F02, and F03.

### Phase 4: scrub and document

- Read every changed file for machine-written tells (banner separators,
  numbered "step" comments, emoji, comments that restate the next line) and
  remove them.
- Confirm no references to AI assistants anywhere in the repository.
- Update `README.md`: drop the stale stages, document the `functions/`
  convention, and list F03 as an active figure.

## Guarantees

- No analysis logic in stages 01 to 03 changes; the seed, path resolution, and
  contrast definitions stay as they are.
- Every deletion is verified unreferenced first; canonical outputs are kept.
- Work lands as small, focused commits, each leaving a clean working tree.

## Out of scope

- Re-running the pipeline end to end. Scripts are edited to preserve behavior;
  a full rerun stays the maintainer's call.
- Renaming `a_script/` to a plural form.
- Changing figure content or statistics.
