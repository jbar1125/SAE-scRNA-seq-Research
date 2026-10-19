#!/usr/bin/env bash
# Continue the declared three-seed CPU experiment, never overwrite result artifacts.
set -euo pipefail
cd /workspace/SAE-scRNA-seq-Research
source /workspace/.venvs/sae-scrna/bin/activate
export XDG_CACHE_HOME=/workspace/.cache MPLCONFIGDIR=/workspace/.cache/matplotlib
export NUMBA_CACHE_DIR=/workspace/.cache/numba
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
run_dir=/workspace/research-runs/norman-sae-primary
active_pid="${1:-}"
if [[ -n "$active_pid" && ! "$active_pid" =~ ^[0-9]+$ ]]; then
  echo 'Expected the PID of the already running seed-0 process' >&2
  exit 2
fi
if [[ "${SAE_DRY_RUN:-0}" != 1 && -n "$active_pid" ]]; then
  while python - "$active_pid" <<'PY'
import pathlib,sys
p=pathlib.Path('/proc')/sys.argv[1]
try:
    args=(p/'cmdline').read_bytes().split(b'\0')
    status=(p/'status').read_text()
    active=any(b'src/program_modulation.py' in a for a in args)
    active=active and any(b'norman-sae-primary/seed0.json' in a for a in args)
    active=active and not any(line.startswith('State:') and 'Z' in line for line in status.splitlines())
except FileNotFoundError:
    active=False
sys.exit(0 if active else 1)
PY
  do
    sleep 30
  done
fi
mkdir -p "$run_dir"
for seed in 0 1 2; do
  result="$run_dir/seed$seed.json"
  checkpoint="$run_dir/seed$seed.pt"
  activations="$run_dir/seed${seed}_activations.npy"
  if [[ -f "$result" ]]; then
    python - "$result" "$seed" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]))
s=r['sae']
assert s['epochs']==s['epochs_completed']==200
assert s['latent_dim']==2048 and s['k']==32 and s['trained_on']=='all'
assert s['n_train_cells']==111445 and s['d']==4081
assert r['provenance']['arguments']['seed']==int(sys.argv[2])
assert r['provenance']['input_sha256']=='a09744916427b08d75dade996d059f717e4812463292a1e5a4e5450589cf2c6a'
assert len(r['programs'])==r['identity_heterogeneity']['n_planned_programs']==9
print('Verified completed seed',sys.argv[2],flush=True)
PY
    continue
  fi
  command=(python -u src/program_modulation.py
    --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad
    --pert-col perturbation --control-value control --train-on all
    --latent 2048 --k 32 --n-hvg 4000 --epochs 200 --batch 1024 --seed "$seed"
    --n-shuffle 0 --n-identity-shuffle 999 --strata-col gemgroup
    --checkpoint-every 5 --checkpoint "$checkpoint"
    --activations-out "$activations" --out "$result")
  if [[ -f "$checkpoint" ]]; then command+=(--resume); fi
  if [[ "${SAE_DRY_RUN:-0}" == 1 ]]; then
    printf '%q ' "${command[@]}"
    printf '\n'
    continue
  fi
  if [[ -e "$activations" ]]; then
    echo "Activation artifact exists without a completed result: $activations; inspect before resuming" >&2
    exit 2
  fi
  echo "Starting/resuming declared full CPU seed $seed" >> "$run_dir/sequence.log"
  "${command[@]}" >> "$run_dir/seed$seed.log" 2>&1
done
