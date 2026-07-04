# CLAUDE.md

Operating manual for this repo. Read top to bottom before doing anything.

**Read [`README.md`](./README.md) for the repo map**, then **[`LAB_NOTEBOOK.md`](./LAB_NOTEBOOK.md)** for the dated day-to-day log (what was done, when, and why). Then the reference docs in [`docs/`](./docs/): **[`PROJECT_HANDOFF.md`](./docs/PROJECT_HANDOFF.md) first** (Phase-1-closeout snapshot), **then [`PROJECT_AUDIT.md`](./docs/PROJECT_AUDIT.md)** (verified inconsistencies, the V0b notebook trap, rigor upgrades, run results — sections J/K are the current findings), **then [`STRATEGY_AND_POSITIONING.md`](./docs/STRATEGY_AND_POSITIONING.md)** (field situating, novelty audit vs Kendiukhov et al., the reframe), **[`PREREGISTRATION.md`](./docs/PREREGISTRATION.md)** (the frozen analysis spec), **[`COMPONENT0_STATUS.md`](./docs/COMPONENT0_STATUS.md)** (Gate-0 state + the runner), and **[`COMPONENT0_RESULTS.md`](./docs/COMPONENT0_RESULTS.md)** (the EXECUTED Gate-0 numbers + frozen SHAs). Where they conflict, the audit wins on facts and the strategy doc wins on framing.

**Repo layout:** `src/` = code, `tests/` = logic tests, `docs/` = narrative/reference docs, `config/` = frozen specs + marker panels (`preregistration_spec.*`, `human_markers.json`), `data/` = raw-input manifest (Drive/large, gitignored), `data_g0/` + `data_g0_human/` = the EXECUTED mouse/human run dirs (large binaries gitignored; small result artifacts committed). Run scripts as `python3 src/<script>.py`; run tests as `python3 tests/<test>.py` (they add `src/` to the path).

---

## 1. CODE DELIVERY RULES (non-negotiable)

**YOU MUST run code before you present it. Never send code you have not executed.** "Looks correct" is not a signal. The only signals are: it compiled, it ran, the test passed.

Before handing over any Python change in this repo, YOU MUST, in order:
1. `python3 -m py_compile <files>` — syntax must pass.
2. `python3 tests/test_v0b_logic.py` — must print `ALL V0b LOGIC TESTS PASSED`.
3. For a script that takes data, do a dry-run or a synthetic-data run that actually executes the code path you changed.
4. Output a verification table and do not claim done until every row says PASS:

   | check | command | result |
   |-------|---------|--------|
   | syntax | py_compile | PASS/FAIL |
   | logic tests | tests/test_v0b_logic.py | PASS/FAIL |
   | executed on real/synthetic data | <command> | PASS/FAIL |

**NEVER** hand the user a command to run before you have run the equivalent yourself. The `torch.manual_seed` NameError shipped on 2026-06-30 happened because this rule was skipped: the refactor removed the top-level `torch` import but `main()` still called it, and it was sent without running. That is the exact failure this section exists to prevent.

If you genuinely cannot run it (missing private data, no GPU), SAY SO explicitly, state what you did verify, and label the rest UNVERIFIED. Do not imply it works.

When code reaches the user's environment (Colab), give a paste that is self-verifying:
- `rm -f <file> <file>.*` then `wget -qO <file> <raw-url>` (plain `wget` writes `<file>.1` and silently runs the stale copy).
- Pin the **commit SHA** in raw GitHub URLs, not the branch, to dodge the ~5-min CDN cache.
- Add a one-line guard that greps for the fix and prints OK/STALE before running.

## 2. HONESTY RULES (non-negotiable)

- **NEVER sugarcoat.** No "great news", no softening a failure into a win. State the result, then the caveat.
- **Distinguish a result from an artifact.** A null caused by an underpowered or blind method is NOT a clean negative. Say which it is. (See PROJECT_AUDIT.md: V0b's per-feature-enrichment test is structurally blind to distributed signal — report that, do not report "erythroid is unified".)
- **Lead with the bad news.** If a finding does not survive a rigorous test, the first sentence says so.
- **Quantify uncertainty.** Provisional, heuristic, and not-yet-statistical numbers must be labeled as such every time.
- **Do not freeze artifacts.** Never SHA-256 "pre-register" an assignment table that came out of a broken or inconclusive run. Diagnose first.
- **Self-correct out loud.** When you ship a bug (you will), name it plainly and fix it. No burying it.
- Prefer brutal constructive criticism over validation. The user's words: "BE AS HONEST AS POSSIBLE, CONSTRUCTIVE CRITICISM TO ALL HEAVENS."

## 3. NON-NEGOTIABLE PROJECT FACTS

- Fully computational (SAEs on scRNA-seq hematopoiesis). No wet bench.
- **Paul15 is pre-log-transformed**: skip `normalize_total`/`log1p`. **CELLxGENE requires both before HVG.** Never conflate the two pipelines.
- Expression matrix MD5 `60183a17983c8b977d036e0f3a58da61` MUST be verified at load (the script asserts it).
- SAEs do NOT preserve feature identity across seeds. Use **per-seed winner-take-all**, NEVER cross-seed feature-index matching.
- Statistical conventions: BH-FDR per test family; Cohen's d with p-values; a finding needs >=3/5 mouse seeds or >=2/3 human seeds.

## 4. PROJECT DEBUG GOTCHAS (learned the hard way)

- **NEVER run the Drive notebooks.** `v0b_module_definitions.ipynb`/`_1`/`_2` — the canonically-named, latest one is v1 with every bug v2 fixed. Use `src/v0b_module_definitions.py` in the repo.
- `torch` is imported lazily inside `load_checkpoints`/`submodule_dynamics`. Do not add a top-level `torch` call; the core logic must run without torch (tests rely on this).
- Inputs live in **My Drive / SAE_scRNA** (consolidated 2026-06-30): `expression_matrix.npy`, `gene_names.csv`, `cell_metadata_palantir.csv`, `sae_seed0-4.pt`. See `data/README.md`.
- The `expression_matrix.npy` is too large to pull through the Drive MCP (base64-into-context). Run where the data already is (Colab), or have the user link-share files for `curl`.

## 5. IMMEDIATE STATUS

- **Gate 0 (Component 0) is EXECUTED end to end and FROZEN** (2026-07-02; see `docs/COMPONENT0_RESULTS.md` for every number and `LAB_NOTEBOOK.md` for the session log). It was run in-container on CPU (SAEs are ~2M params on 2-4k cells; no GPU needed). Both arms done:
  - Mouse (Paul15, 2730x2012): 10-seed 512-latent SAE, l1=0.8, L0 mean 35.75 (in 20-50 band). Original asymmetric-modularity claim **NOT SUPPORTED** (0/10 seeds); the data show the REVERSE (granulocyte 3 real programs, erythroid 1 dominant hemoglobin program) and it is seed-robust. Frozen bundle sha256 `172861...`.
  - Human (Setty 2019 CD34+ marrow, 4142x2032): 10-seed SAE, L0 mean 32.5. Claim NOT SUPPORTED (0/10); granulocyte-more-distributed direction replicates; human globin program not significant vs the abundance-matched null (CD34+ progenitor selection). Frozen bundle sha256 `cc8159c8...`.
  - **The key finding is method-dependence:** PCA/NMF/SAE give different modularity answers on the identical matrix; the original claim survives in exactly 1 of 6 method x species combinations (mouse NMF). Marker leave-one-out: 0/34 (mouse), 0/36 (human) flips. Null calibrated (Progenitor non-significant both species).
- The prior "erythroid undetected" runs (Audit F-I) were a marker-curation artifact (wrong globin symbols, missing Alas2), fixed. Do not cite them as the finding.
- Reframe (STRATEGY_AND_POSITIONING.md): the novel contribution is the artifact-vs-signal framework + the causal CRISPRi test, NOT "SAEs on single-cell" (now published: Kendiukhov 2026 arXiv 2603.02952 and the 2025-2026 wave). Lead with the framework and the causal spine.
- **Next: Component 2** (Replogle CRISPRi causal head-to-head vs the 6.2% embedding-space null) is the centerpiece; optional cheap hardening first (SCENIC as a 5th baseline, TRRUST/ChIP-Atlas cross-reference). Plan in `docs/STRATEGY_AND_POSITIONING.md` and `docs/phase2_research_plan_v6.md`. Do NOT re-freeze the frozen Gate-0 bundles; new work is new components.

## 6. WORKING STYLE

Extremely concise. No em-dashes. No filler. Per-bullet formatting. Code ready-to-run. First person fine in operational contexts. ROL voice: formal third-person. Methodology voice: imperative future. See handoff section 15.

---

## SMART SUMMARY (read this if you read nothing else)

This project's job is to find out whether SAE-derived gene programs in hematopoiesis are real and rigorous, not to confirm a story. The single most important behavior: **run the code, prove it works, then report the result without spin.** The headline claim ("asymmetric modularity": granulocyte unified, erythroid distributed) currently has **no valid quantitative support** — the v1 notebook that "supported" it used high-abundance markers v2 correctly drops, and the rigorous v2 test came back inconclusive because it cannot detect distributed signal at all. Treat every prior number as provisional until re-derived in version-controlled, executed, tested code. When you write code, debug it before sending; when you report findings, lead with what failed. If you are ever choosing between sounding helpful and being accurate, choose accurate.
