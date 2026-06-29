# CLAUDE.md

**Read [`PROJECT_HANDOFF.md`](./PROJECT_HANDOFF.md) first.** It is the single source of truth for project state at Phase 1 closeout (15 sections: thesis, findings, V0-V15 plan, active bugs, Phase 2 plan, file inventory, references, next actions).

## Non-negotiable facts
- Project is fully computational (SAEs on scRNA-seq hematopoiesis). No wet bench.
- **Paul15 is pre-log-transformed**: skip `normalize_total`/`log1p`. **CELLxGENE requires both before HVG.** Do not conflate the two pipelines.
- Expression matrix MD5: `60183a17983c8b977d036e0f3a58da61` (verify at the start of every notebook).
- SAEs do NOT preserve feature identity across seeds. Use per-seed winner-take-all, never cross-seed feature-index matching.
- Pre-register thresholds via SHA-256 freezes before inspecting data.

## Immediate blocker
Run V0b v2 (`v0b_module_definitions.ipynb`, verify v2 not v1 per handoff 4.7) to completion and freeze `module_assignments_v0b.csv`. Every downstream Phase 2 claim depends on whether asymmetric modularity survives the gene-content test.

## Working style
Extremely concise. No em-dashes. No filler. Brutal honest methodological self-assessment over validation. Code delivered ready-to-run. See handoff section 15.
