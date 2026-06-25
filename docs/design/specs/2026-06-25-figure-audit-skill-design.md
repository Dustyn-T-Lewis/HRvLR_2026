# Figure Audit Skill — Design Spec

**Date:** 2026-06-25
**Status:** Approved design, pending implementation plan
**Artifact:** `.claude/skills/figure-audit/SKILL.md` (project-scoped reusable skill)

## Purpose

A reusable skill that reviews the HRvLR figures and panels against the study
hypothesis — does each figure communicate and answer the questions it is meant
to — and then enforces a clean, consistent directory and data layout. It is a
thinking partner: it interrogates intent, grounds the review in the literature,
reviews scope, and only then reorganizes, with a gate before any change.

Invocation: `/figure-audit` (all figures) or `/figure-audit F02` (one figure).

## Anchor: hypothesis and figure roles (encoded in the skill)

**Hypothesis.** In resistance-trained skeletal muscle, high vs low hypertrophy
responders diverge in their proteome — characterized phenotypically (F01), in
the differential and enriched proteome (F02), and mechanistically and
predictively via co-expression modules linked to the adaptation (F03).

| Figure | Role | Question it must answer |
| --- | --- | --- |
| F01 | Descriptive overview | Who are the responders? HR vs LR phenotype (growth/strength) despite matched training. |
| F02 | Proteome overview + enrichment volcanoes | What is different? DEPs per contrast plus enriched biology (enrichVolcano rings). |
| F03 | Clustering — mechanistic + predictive | How, and can we predict? Modules linked to adaptation (mechanism) and baseline proteome predicting the response (prediction). |

The role table is the skill's reference; it is extended if figures are added.

## What the skill produces

1. **Figure-level narrative check** — does each figure answer its question, and
   do F01 to F03 build one coherent argument with no gaps or redundancy?
2. **Per-panel verdict table** — one row per panel:
   `panel | what it shows | hypothesis-question served | role-fit | verdict | action`,
   where verdict is one of KEEP, CUT, REFRAME, MERGE, ADD-MISSING.

## Organization standard the skill enforces (per figure)

- `a_script/` — ordered setup/run/panel scripts; a `supp/` subdirectory for
  supplementary-only scripts.
- `b_reports/` — the composite figure at the top level, a `panels/`
  subdirectory for individual panels, and a `supp/` subdirectory for
  supplementary panels.
- `c_data/` — one concise, well-named Excel workbook that is both the figure's
  results (output) and a clean input for downstream scripts and panels: one
  sheet per panel's underlying data plus a metadata sheet, tidy columns, no
  orphan or duplicate sheets. Only the RDS/CSV genuinely consumed downstream are
  kept alongside it.

## Workflow: focus, research, review, apply

**Phase 0 — Focusing questions (brainstorm first).** Before judging anything,
ask targeted questions one at a time: the single claim each figure must land;
the target venue/audience (sets the bar); what is explicitly out of scope; any
panel already suspected weak or redundant. Propose the framing back and get a
nod before reviewing.

**Phase 1 — Literature review.** For each figure, check against published work
using the research-agent and MCP pattern (PubMed, OpenAlex, Semantic Scholar;
verified citations): are the analyses canonical and defensible; what do
comparable responder/hypertrophy-proteomics papers show in their figures; what
will a reviewer expect, and what are the known pitfalls. Flag claims the data
cannot support and gaps the literature says should be covered.

**Phase 2 — Scope review.** Produce the figure-level narrative and the per-panel
verdict table, informed by Phases 0 and 1.

**Phase 3 — Gated apply.** Present the plan and STOP for approval. On approval:
reframe, cut, or merge panels; move files into the standard layout; build or
refresh the Excel; re-run scripts and verify (panels render; `c_data` byte
unchanged where it should not move).

## Guardrails

- Snapshot before any move or delete; everything stays git-recoverable.
- After reorganization and re-run, verify renders and data integrity.
- Nothing changes before the user's sign-off (Phase 3 gate).

## Out of scope

- Generating new analyses or panels beyond ADD-MISSING recommendations the user
  approves at the gate.
- Editing manuscript prose; the skill audits figures, not the paper text.
