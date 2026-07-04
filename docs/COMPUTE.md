# COMPUTE.md — what runs where, and why

A standing record of the compute decisions, so a mentor or judge can see the
hardware reasoning and reproduce the runs. Short version: **Gate 0 needs no GPU and
was run on CPU; the GPU requirement starts at Component 2.**

## The environments in play

| Environment | What it is | GPU |
|-------------|-----------|-----|
| This dev container | Claude Code on the web: an ephemeral, isolated Linux cloud VM, internet via a proxy. Where Gate 0 was executed. | **None** (`torch.cuda.is_available()` is False; 4 CPU cores) |
| Colab (T4) | NVIDIA T4, CUDA-capable. The intended GPU path for Component 2. | NVIDIA T4 |
| A local PC with RX 5700 XT | AMD RDNA1 (Navi 10 / gfx1010). | AMD (see caveat) |

## Why Gate 0 ran on CPU (and that is fine)

- The SAE is ~2.06M parameters: `Linear(2012->512)` encoder + `Linear(512->2012)`
  decoder. Tiny.
- The data is small: mouse 2730 x 2012, human 4142 x 2032 (~20-34 MB, fits in RAM).
- Measured timing on 4 CPU cores: ~50 s/seed (mouse), ~75 s/seed (human) at 300
  epochs. A full 10-seed run is ~8-12 min. L1 sweeps are ~50 s each.
- CPU and GPU compute the same math; with fixed seeds any difference is ~1e-6
  floating-point rounding and changes no verdict. So CPU is a correctness-neutral,
  speed-only choice here, not a compromise.
- Honest caveat: CPU-only means the CUDA code path itself was not exercised in this
  environment. The numeric result is identical because the model/seeds are fixed.

## Why the RX 5700 XT is not a practical GPU for this project

- It is an **AMD** card. PyTorch `torch.cuda` acceleration is **NVIDIA/CUDA only**.
  AMD needs the **ROCm** backend instead.
- The 5700 XT is RDNA1 / gfx1010, which is **not officially supported** by ROCm for
  PyTorch. It sometimes works via `HSA_OVERRIDE_GFX_VERSION` hacks, but it is
  fragile and frequently does not work at all.
- Conclusion: do not plan on the 5700 XT for the ML stack. For the GPU stage use an
  **NVIDIA** card or the Colab **T4** (CUDA), which works out of the box.

## What actually needs a GPU (starts at Component 2)

- **Component 2** (the causal centerpiece): Replogle CRISPRi head-to-head vs
  Kendiukhov's 6.2% embedding-space null. Involves large single-cell foundation-model
  embeddings (scGPT-scale). Needs real GPU memory. NOT runnable in this CPU container.
- A **large human atlas** (25k+ CELLxGENE cells) instead of the 4142-cell Setty set.

## How to run the GPU stage

- **Colab (simplest):** T4 runtime, `pip install` the deps, clone the repo, run the
  Component 2 scripts. CUDA is automatic.
- **Own NVIDIA hardware via Remote Control:** run `claude remote-control` (or
  `--remote-control`) locally so Claude Code executes on your machine with your GPU,
  filesystem, and data staying local, and steer/monitor from claude.ai or the phone
  app. See `docs/` / the Remote Control docs.

## Operational risk that is real (and unrelated to CPU vs GPU)

This dev container is **ephemeral**: it has restarted mid-session and killed
in-progress background jobs. Mitigation, always applied: commit + push after every
completed unit, and keep large regenerable binaries out of git. Nothing that matters
lives only in the container.
