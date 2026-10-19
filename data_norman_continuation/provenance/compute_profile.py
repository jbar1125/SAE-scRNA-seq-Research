#!/usr/bin/env python3
"""Compute-only Norman primary-shape timing; no expression outcomes are reported."""
import os
for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

from pathlib import Path
import hashlib
import json
import math
import platform
import time

import anndata as ad
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

path = Path('/workspace/datasets/sae-scrna/norman_raw.h5ad')
out = Path('/workspace/research-drafts/norman_preflight/compute_profile.json')
panel_width = 4081
batch_size = 1024
latent = 2048
top_k = 32
epochs = 200
warmup = 3
repeats = 8
a = ad.read_h5ad(path, backed='r')
rows = a.X[:batch_size, :panel_width]
rep = np.asarray(rows.toarray() if hasattr(rows, 'toarray') else rows, dtype=np.float32)
n_cells = a.n_obs
a.file.close()
torch.set_num_interop_threads(1)
measurements = []
for threads in (1, 2, 4):
    torch.set_num_threads(threads)
    torch.manual_seed(0)
    enc = nn.Linear(panel_width, latent)
    dec = nn.Linear(latent, panel_width, bias=False)
    with torch.no_grad():
        dec.weight.copy_(F.normalize(dec.weight, dim=0))
    opt = torch.optim.AdamW(list(enc.parameters()) + list(dec.parameters()), lr=4e-4)
    matrix = torch.tensor(rep - rep.mean(0, keepdims=True), dtype=torch.float32)

    def step():
        perm = torch.randperm(matrix.shape[0])
        xb = matrix[perm]
        values, index = torch.topk(enc(xb), top_k, dim=1)
        activation = torch.zeros((len(xb), latent), dtype=xb.dtype)
        activation.scatter_(1, index, F.relu(values))
        loss = F.mse_loss(dec(activation), xb)
        opt.zero_grad()
        loss.backward()
        opt.step()
        with torch.no_grad():
            dec.weight.copy_(F.normalize(dec.weight, dim=0))

    for _ in range(warmup):
        step()
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        step()
        samples.append(time.perf_counter() - started)
    median = float(np.median(samples))
    measurements.append({
        'threads': threads, 'step_seconds': samples,
        'median_seconds': median, 'min_seconds': min(samples), 'max_seconds': max(samples),
        'all_cell_training_hours_estimate': median * math.ceil(n_cells / batch_size) * epochs / 3600,
        'control_only_training_hours_estimate': median * math.ceil(11855 / batch_size) * epochs / 3600,
    })
    print(json.dumps(measurements[-1]), flush=True)
    del opt, enc, dec, matrix

result = {
    'scope': 'Timing only. No loss, expression differences, or biological results were inspected or reported.',
    'input': str(path), 'input_bytes': path.stat().st_size,
    'profile_rows': batch_size, 'profile_columns': panel_width,
    'column_selection': 'First4081 raw genes, not data-selected HVGs; used only to profile full matrix shape',
    'caveat': 'Input is a small raw-count slice. Full workflow log-normalizes and chooses HVGs; arithmetic timing should be comparable, not a validated training run.',
    'configuration': {'latent': latent, 'k': top_k, 'batch': batch_size, 'epochs': epochs,
                      'warmup_steps': warmup, 'timed_steps': repeats, 'optimizer': 'AdamW', 'lr': 4e-4},
    'n_all_cells': n_cells, 'all_cell_steps': math.ceil(n_cells / batch_size) * epochs,
    'control_only_steps': math.ceil(11855 / batch_size) * epochs,
    'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'torch': torch.__version__},
    'affinity_cpus': len(os.sched_getaffinity(0)),
    'cgroup_cpu_max': Path('/sys/fs/cgroup/cpu.max').read_text().strip(),
    'cgroup_memory_max': Path('/sys/fs/cgroup/memory.max').read_text().strip(),
    'profile_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'measurements': measurements,
}
out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'saved': str(out), 'measurements': measurements}), flush=True)
