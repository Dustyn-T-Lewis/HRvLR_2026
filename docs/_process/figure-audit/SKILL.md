---
name: figure-audit
description: Review the HRvLR manuscript figures and panels against the study hypothesis (figure-level narrative + per-panel verdicts), ground the review in the literature, then enforce a clean a_script/b_reports/c_data layout. Use when reviewing figure scope, deciding keep/cut/reframe, or reorganizing a figure directory. Runs focus -> research -> review -> apply, gated before any change. Invoke /figure-audit or /figure-audit F02.
---

# Figure Audit

Reviews the HRvLR figures against the hypothesis and reorganizes them to a clean
standard. A thinking partner, not a linter: it asks focusing questions, grounds
the review in the literature, reviews scope, then applies — with a gate before
any file changes. Invoke `/figure-audit` (all figures) or `/figure-audit F02`
(one figure).

## Anchor: hypothesis and figure roles

**Hypothesis.** In resistance-trained skeletal muscle, high vs low hypertrophy
responders diverge in the magnitude and trajectory of their proteomic training
response, not merely in baseline state — characterized phenotypically (F01), in
the differential and enriched proteome (F02, F03), in HR-vs-LR concordance of the
response (F04, F05), and via co-expression modules linked to the adaptation (F06).

| Figure | Role | Question it must answer |
| --- | --- | --- |
| F01 | Descriptive overview | Who are the responders? HR vs LR phenotype (growth/strength) despite matched training. |
| F02 | Proteome overview + QC | What is different? DEPs per contrast, effect sizes, variance structure. |
| F03 | enrichVolcano ring-volcanoes | What moves where? logFC/pi volcanoes with fgsea NES rings per contrast; also builds the shared fgsea cache. |
| F04 | Training-phase concordance | Do HR and LR move the same biology the same way from T1 to T2, or diverge? |
| F05 | Acute-phase concordance | Same question for the acute bout, T2 to T3. |
| F06 | WGCNA module-phenotype | Which co-expression modules track the adaptation (within-cohort association)? |

Extend this table if figures are added.

## Organization standard (enforced per figure)

- `a_script/` — ordered setup/run/panel scripts; a `supp/` subdirectory for
  supplementary-only scripts.
- `b_reports/` — the composite figure at the top level, a `panels/`
  subdirectory for individual panels, a `supp/` subdirectory for supplementary
  panels.
- `c_data/` — one concise, well-named Excel workbook that is both the figure's
  results (output) and a clean input for downstream scripts/panels: one sheet
  per panel's underlying data plus a metadata sheet, tidy columns, no orphan or
  duplicate sheets. Keep only the RDS/CSV genuinely consumed downstream.

## Workflow: focus -> research -> review -> apply

### Phase 0 — Focusing questions (ask first, one at a time)

Before judging anything, ask one question per message:
- What single claim must each target figure land?
- Target venue/audience (sets the bar)?
- What is explicitly out of scope (so we cut, not just critique)?
- Any panel already suspected to be weak or redundant?

Propose the framing back and get a nod before reviewing.

### Phase 1 — Literature review

For each target figure, dispatch research — the `deep-research` skill, or the
project's literature MCP servers (PubMed, OpenAlex, Semantic Scholar) via
research subagents; verified citations only — flag anything unverifiable:
- Methods: are the analyses canonical and defensible?
- Framing and expectations: what do comparable responder / hypertrophy-proteomics
  papers show in their figures; what will a reviewer expect; known pitfalls.

Flag claims the data cannot support and gaps the literature says to cover.

### Phase 2 — Scope review (produce both)

1. **Figure-level narrative check** — does each figure answer its assigned
   question, and do F01 -> F02 -> F03 build one coherent argument with no gaps
   or redundancy?
2. **Per-panel verdict table** — one row per panel:

   | panel | what it shows | hypothesis-question served | role-fit | verdict | action |

   `verdict` is one of KEEP, CUT, REFRAME, MERGE, ADD-MISSING; `action` is the
   one-line change implied by the verdict.

### Phase 3 — Gated apply

Present the narrative + verdict table + reorganization plan, then STOP for
approval. On approval:
- Reframe, cut, or merge panels per the verdicts.
- Move files into the standard layout (a_script/supp, b_reports composite +
  panels + supp).
- Build or refresh the `c_data` Excel.
- Re-run scripts and verify: panels render, and `c_data` is byte unchanged where
  it should not move.

## Guardrails

- Snapshot before any move or delete; everything stays git-recoverable.
- After reorganization and re-run, verify renders and data integrity (byte-diff
  the `c_data` snapshots; visually check a render — never md5 an image).
- Nothing changes before the user's sign-off at the Phase 3 gate.
- `/figure-audit F0x` audits only that figure; with no argument, audit all and
  also check the F01 -> F02 -> F03 narrative arc.

## Out of scope

- Generating new analyses or panels beyond ADD-MISSING recommendations the user
  approves at the gate.
- Editing manuscript prose; the skill audits figures, not the paper text.
