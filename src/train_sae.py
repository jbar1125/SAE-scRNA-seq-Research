#!/usr/bin/env python3
"""Train sparse autoencoders on the (scaled) expression matrix.

Compute/data upgrade: supports OVERCOMPLETE dictionaries (latent_dim >> 128) so
distributed programs have room to separate into distinct features. The 128-latent
compressive default is kept as a baseline. Checkpoints are saved as plain
state_dicts, which v0b_module_definitions.load_checkpoints reads (it infers
latent_dim from the decoder shape, so overcomplete loads with no other change).

Design follows PROJECT_HANDOFF.md 4.3: Linear(D->L) encoder, ReLU, Linear(L->D)
decoder; AdamW lr 5e-4; L1 penalty lambda=0.1; decoder columns unit-normalized
each step; dead-latent resampling; 300 epochs; batch 256. Reports recon MSE and L0
(mean active latents/cell) per seed; the 20-50 L0 band is the project QC gate.

Run
  python train_sae.py --matrix expression_matrix.npy --latent-dim 512 \
      --out-dir checkpoints_L512 --seeds 5
Produces sae_seed0.pt ... sae_seed{S-1}.pt + a training_report.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class SparseAutoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super().__init__()
        self.encoder = nn.Linear(input_dim, latent_dim, bias=True)
        self.decoder = nn.Linear(latent_dim, input_dim, bias=True)

    def forward(self, x):
        z = F.relu(self.encoder(x))
        return self.decoder(z), z

    @torch.no_grad()
    def normalize_decoder(self):
        # unit-normalize decoder columns (per-latent dictionary atoms)
        w = self.decoder.weight  # (input_dim, latent_dim)
        self.decoder.weight.copy_(w / (w.norm(dim=0, keepdim=True) + 1e-8))


def train_one(X, latent_dim, seed, epochs, batch, lr, l1, device,
              resample_every=500, dead_frac_window=0.0):
    torch.manual_seed(seed)
    np.random.seed(seed)
    n, d = X.shape
    model = SparseAutoencoder(d, latent_dim).to(device)
    model.normalize_decoder()
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-5)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    Xt = torch.tensor(X, dtype=torch.float32, device=device)

    step = 0
    fired = torch.zeros(latent_dim, device=device)  # activation counts since last resample
    for _ in range(epochs):
        perm = torch.randperm(n, device=device)
        for i in range(0, n, batch):
            xb = Xt[perm[i:i + batch]]
            recon, z = model(xb)
            loss = F.mse_loss(recon, xb) + l1 * z.abs().mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            model.normalize_decoder()
            fired += (z > 0).float().sum(0)
            step += 1
            if step % resample_every == 0:
                dead = (fired == 0)
                if dead.any():
                    # re-init dead latents' encoder rows toward high-error directions
                    with torch.no_grad():
                        err = (recon - xb)
                        # pick random high-error samples as new directions
                        idx = torch.randint(0, xb.shape[0], (int(dead.sum()),), device=device)
                        newdir = F.normalize(err[idx], dim=1)
                        model.encoder.weight[dead] = newdir * 0.1
                        model.encoder.bias[dead] = 0.0
                fired.zero_()
        sched.step()

    model.eval()
    with torch.no_grad():
        recon, z = model(Xt)
        mse = F.mse_loss(recon, Xt).item()
        l0 = (z > 0).float().sum(1).mean().item()
    return model, {"recon_mse": mse, "mean_active_l0": l0}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True, help="path to scaled expression_matrix .npy")
    ap.add_argument("--out-dir", default="./checkpoints")
    ap.add_argument("--latent-dim", type=int, default=128)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--lr", type=float, default=5e-4)
    ap.add_argument("--l1", type=float, default=0.1)
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    X = np.load(args.matrix).astype(np.float32)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    print(f"X {X.shape} | latent_dim {args.latent_dim} | device {device}")

    report = {"latent_dim": args.latent_dim, "input_dim": X.shape[1],
              "n_cells": X.shape[0], "epochs": args.epochs, "seeds": {}}
    for seed in range(args.seeds):
        model, stats = train_one(X, args.latent_dim, seed, args.epochs,
                                 args.batch, args.lr, args.l1, device)
        torch.save(model.state_dict(), out / f"sae_seed{seed}.pt")
        report["seeds"][seed] = stats
        flag = "OK" if 20 <= stats["mean_active_l0"] <= 50 else "OUTSIDE 20-50 QC BAND"
        print(f"seed {seed}: recon_mse {stats['recon_mse']:.4f} | "
              f"L0 {stats['mean_active_l0']:.1f} [{flag}]")
    with open(out / "training_report.json", "w") as f:
        json.dump(report, f, indent=2)
    print(f"Saved {args.seeds} checkpoints to {out}")


if __name__ == "__main__":
    main()
