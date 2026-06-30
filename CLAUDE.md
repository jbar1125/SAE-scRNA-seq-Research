# CLAUDE.md

**Read [`PROJECT_HANDOFF.md`](./PROJECT_HANDOFF.md) first** (Phase-1-closeout snapshot: thesis, findings, V0-V15 plan, active bugs, Phase 2 plan, file inventory, references, next actions), **then [`PROJECT_AUDIT.md`](./PROJECT_AUDIT.md)** for the correction + rigor layer (verified inconsistencies, the V0b notebook trap, rigor upgrades). Where they conflict, the audit wins.

## Non-negotiable facts
- Project is fully computational (SAEs on scRNA-seq hematopoiesis). No wet bench.
- **Paul15 is pre-log-transformed**: skip `normalize_total`/`log1p`. **CELLxGENE requires both before HVG.** Do not conflate the two pipelines.
- Expression matrix MD5: `60183a17983c8b977d036e0f3a58da61` (verify at the start of every notebook).
- SAEs do NOT preserve feature identity across seeds. Use per-seed winner-take-all, never cross-seed feature-index matching.
- Pre-register thresholds via SHA-256 freezes before inspecting data.

## Immediate blocker
Run the canonical v2 V0b: `python v0b_module_definitions.py --data-dir <inputs>` (NOT the Drive notebooks; the latest-named one is v1, see PROJECT_AUDIT.md B). Then commit the frozen `module_assignments_v0b.csv` + `v0b_provenance.json`. Every downstream Phase 2 claim depends on whether asymmetric modularity survives the gene-content test, and it has not been run to a frozen result yet.

## Working style
Extremely concise. No em-dashes. No filler. Brutal honest methodological self-assessment over validation. Code delivered ready-to-run. See handoff section 15.
