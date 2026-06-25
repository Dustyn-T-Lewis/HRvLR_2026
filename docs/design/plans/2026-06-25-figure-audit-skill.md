# Figure Audit Skill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a project-scoped reusable skill, `/figure-audit`, that reviews the HRvLR figures against the study hypothesis (figure-level narrative + per-panel verdicts), grounds the review in the literature, and enforces a clean `a_script`/`b_reports`/`c_data` layout — gated before any change.

**Architecture:** A single self-contained `SKILL.md` at `.claude/skills/figure-audit/` in the repo. YAML frontmatter makes it discoverable as `/figure-audit`; the body encodes the hypothesis, the figure-role table, the focus→research→review→apply workflow, the directory/Excel standard, and guardrails. No code, no runtime — the artifact is an instruction document the model follows.

**Tech Stack:** Markdown + YAML frontmatter (Claude Code skill format). Verification via `python3 -c` YAML parse and `grep` section checks.

## Global Constraints

- Artifact path: `.claude/skills/figure-audit/SKILL.md` (project-scoped, travels with the repo).
- Invocation: `/figure-audit` (all figures) or `/figure-audit F02` (one figure).
- Spec source of truth: `docs/design/specs/2026-06-25-figure-audit-skill-design.md`.
- Verdict categories (exact): KEEP, CUT, REFRAME, MERGE, ADD-MISSING.
- Standard layout: `a_script/` (+`supp/`), `b_reports/` (composite + `panels/` + `supp/`), `c_data/` (one concise Excel = output + clean input).
- Workflow order (exact): focus → research → review → apply, gated before any change.
- Commit style: one short sentence, no AI/attribution trailers.

---

### Task 1: Scaffold skill directory and frontmatter

**Files:**
- Create: `.claude/skills/figure-audit/SKILL.md`

**Interfaces:**
- Produces: a discoverable skill named `figure-audit` with valid YAML frontmatter (`name`, `description`).

- [ ] **Step 1: Create the file with frontmatter only**

Create `.claude/skills/figure-audit/SKILL.md` with exactly:

```markdown
---
name: figure-audit
description: Review the HRvLR manuscript figures and panels against the study hypothesis (figure-level narrative + per-panel verdicts), ground the review in the literature, then enforce a clean a_script/b_reports/c_data layout. Use when reviewing figure scope, deciding keep/cut/reframe, or reorganizing a figure directory. Runs focus -> research -> review -> apply, gated before any change. Invoke /figure-audit or /figure-audit F02.
---

# Figure Audit
```

- [ ] **Step 2: Verify the frontmatter parses as YAML**

Run:
```bash
python3 -c "import yaml,sys; d=yaml.safe_load(open('.claude/skills/figure-audit/SKILL.md').read().split('---')[1]); assert d['name']=='figure-audit'; assert len(d['description'])>40; print('frontmatter OK')"
```
Expected: `frontmatter OK`

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/figure-audit/SKILL.md
git commit -m "scaffold the figure-audit skill"
```

---

### Task 2: Write the anchor and organization standard

**Files:**
- Modify: `.claude/skills/figure-audit/SKILL.md` (append body after `# Figure Audit`)

**Interfaces:**
- Consumes: the frontmatter + `# Figure Audit` heading from Task 1.
- Produces: the `## Anchor` and `## Organization standard` sections later tasks reference.

- [ ] **Step 1: Append the intro, anchor, and organization sections**

Append to `.claude/skills/figure-audit/SKILL.md`:

```markdown

Reviews the HRvLR figures against the hypothesis and reorganizes them to a clean
standard. A thinking partner, not a linter: it asks focusing questions, grounds
the review in the literature, reviews scope, then applies — with a gate before
any file changes. Invoke `/figure-audit` (all figures) or `/figure-audit F02`
(one figure).

## Anchor: hypothesis and figure roles

**Hypothesis.** In resistance-trained skeletal muscle, high vs low hypertrophy
responders diverge in their proteome — characterized phenotypically (F01), in
the differential and enriched proteome (F02), and mechanistically and
predictively via co-expression modules linked to the adaptation (F03).

| Figure | Role | Question it must answer |
| --- | --- | --- |
| F01 | Descriptive overview | Who are the responders? HR vs LR phenotype (growth/strength) despite matched training. |
| F02 | Proteome overview + enrichment volcanoes | What is different? DEPs per contrast plus enriched biology (enrichVolcano rings). |
| F03 | Clustering — mechanistic + predictive | How, and can we predict? Modules linked to adaptation (mechanism) and baseline proteome predicting the response (prediction). |

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
```

- [ ] **Step 2: Verify the anchor and standard sections exist**

Run:
```bash
grep -qc "## Anchor: hypothesis and figure roles" .claude/skills/figure-audit/SKILL.md && grep -q "## Organization standard" .claude/skills/figure-audit/SKILL.md && grep -q "enrichVolcano rings" .claude/skills/figure-audit/SKILL.md && echo "anchor+standard OK"
```
Expected: `anchor+standard OK`

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/figure-audit/SKILL.md
git commit -m "add the figure-audit anchor and layout standard"
```

---

### Task 3: Write the workflow, output format, and guardrails

**Files:**
- Modify: `.claude/skills/figure-audit/SKILL.md` (append after the Organization section)

**Interfaces:**
- Consumes: the `## Anchor` (figure roles) and `## Organization standard` from Task 2.
- Produces: the complete, runnable skill body (workflow + verdict table + guardrails).

- [ ] **Step 1: Append the workflow, output format, and guardrails**

Append to `.claude/skills/figure-audit/SKILL.md`:

```markdown

## Workflow: focus -> research -> review -> apply

### Phase 0 — Focusing questions (ask first, one at a time)

Before judging anything, ask one question per message:
- What single claim must each target figure land?
- Target venue/audience (sets the bar)?
- What is explicitly out of scope (so we cut, not just critique)?
- Any panel already suspected to be weak or redundant?

Propose the framing back and get a nod before reviewing.

### Phase 1 — Literature review

For each target figure, dispatch research subagents / MCP tools (PubMed,
OpenAlex, Semantic Scholar; verified citations only — flag anything unverifiable):
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
```

- [ ] **Step 2: Verify all workflow phases and verdict categories are present**

Run:
```bash
f=.claude/skills/figure-audit/SKILL.md
grep -q "### Phase 0" $f && grep -q "### Phase 1" $f && grep -q "### Phase 2" $f && grep -q "### Phase 3" $f && grep -q "KEEP, CUT, REFRAME, MERGE, ADD-MISSING" $f && grep -q "STOP for" $f && echo "workflow OK"
```
Expected: `workflow OK`

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/figure-audit/SKILL.md
git commit -m "add the figure-audit workflow, verdicts, and guardrails"
```

---

### Task 4: Validate the complete skill against the spec

**Files:**
- Read: `.claude/skills/figure-audit/SKILL.md`, `docs/design/specs/2026-06-25-figure-audit-skill-design.md`

**Interfaces:**
- Consumes: the complete SKILL.md from Tasks 1–3.
- Produces: a validated, spec-complete skill (no code interface).

- [ ] **Step 1: Confirm every spec requirement maps to a SKILL.md section**

Run:
```bash
f=.claude/skills/figure-audit/SKILL.md
for s in "Hypothesis" "Descriptive overview" "enrichment volcanoes" "mechanistic" "predictive" "Focusing questions" "Literature review" "verdict table" "Organization standard" "Guardrails" "Gated apply"; do
  grep -qi "$s" "$f" || echo "MISSING: $s"
done
echo "spec-coverage check done"
```
Expected: `spec-coverage check done` with no `MISSING:` lines.

- [ ] **Step 2: Confirm the skill is discoverable (frontmatter name + body present)**

Run:
```bash
head -5 .claude/skills/figure-audit/SKILL.md | grep -q "name: figure-audit" && test $(wc -l < .claude/skills/figure-audit/SKILL.md) -gt 60 && echo "skill complete"
```
Expected: `skill complete`

- [ ] **Step 3: Commit any fixes**

```bash
git add .claude/skills/figure-audit/SKILL.md
git commit -m "validate the figure-audit skill against the spec" --allow-empty
```

---

## Self-Review

- **Spec coverage:** anchor/hypothesis (Task 2), figure roles (Task 2), both outputs — narrative + per-panel verdicts (Task 3), organization standard incl. Excel (Task 2), focus→research→review→apply workflow (Task 3), gated apply + guardrails (Task 3), per-figure scope (Task 3). Task 4 grep-verifies all spec keywords map to a section. No gaps.
- **Placeholder scan:** every task contains the literal markdown to write and exact verify commands; no TBD/TODO.
- **Type consistency:** verdict categories `KEEP, CUT, REFRAME, MERGE, ADD-MISSING` and the layout paths are identical in the Global Constraints, Task 2, and Task 3.
