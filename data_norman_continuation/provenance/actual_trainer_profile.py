#!/usr/bin/env python3
"""Shape-only timing of the actual trainer; no biological observations."""
import os
for name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
import hashlib
import json
import math
import sys
import time
import numpy as np
import torch
sys.path.insert(0,'/workspace/SAE-scRNA-seq-Research/src')
import causal_pipeline as cp

torch.set_num_interop_threads(1)
torch.set_num_threads(1)
rng=np.random.default_rng(0)
# Warm optimizer/module initialization outside all reported timings.
encode,_=cp.train_topk_sae(rng.normal(size=(128,64)).astype(np.float32),32,8,epochs=1,batch=128,device='cpu')
del encode
matrix=rng.normal(size=(4096,4081)).astype(np.float32)
rows=[]
for threads in (2,4):
    torch.set_num_threads(threads)
    for latent in (2048,512):
        values={}
        for epochs in (2,4):
            started=time.perf_counter()
            encode,info=cp.train_topk_sae(matrix,latent,32,epochs=epochs,batch=1024,device='cpu')
            elapsed=time.perf_counter()-started
            del encode
            values[epochs]=elapsed
        slope=(values[4]-values[2])/8
        coarse=values[2]/8
        point={'threads':threads,'latent':latent,'fit2epochs8steps_seconds':values[2],
               'fit4epochs16steps_seconds':values[4],
               'incremental_seconds_per_step':slope,
               'coarse_seconds_per_step':coarse,
               'primary_200epochs_hours':slope*21800/3600 if latent==2048 else None,
               'pilot_50epochs_minutes':slope*5450/60 if latent==512 else None,
               'three_pilots_50epochs_minutes':slope*5450*3/60 if latent==512 else None}
        rows.append(point)
        print(json.dumps(point),flush=True)
output={'scope':'Actual trainer on synthetic random inputs; timing only, no biological analysis',
        'shape':[4096,4081],'k':32,'batch':1024,'projection_all_cells':111445,
        'projection_steps_per_epoch':109,'models':rows,
        'caveat':'Short fits and differences estimate arithmetic only; scheduler contention, loading, HVG, encoding, checkpoint IO, and permutations are excluded. Not a runtime guarantee.',
        'torch':torch.__version__,'trainer_source_sha256':hashlib.sha256(Path(cp.__file__).read_bytes()).hexdigest(),
        'profiler_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path('/workspace/research-drafts/norman_preflight/actual_trainer_profile.json').write_text(json.dumps(output,indent=2)+'\n')
