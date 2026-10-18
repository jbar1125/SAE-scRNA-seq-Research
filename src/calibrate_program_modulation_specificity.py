#!/usr/bin/env python3
"""Calibrate the conditional identity test on independent synthetic null datasets.

Default, fixed audit: 200 independent null simulations, nine distinct matched
features/programs, 100 controls, ten single-perturbation groups of 40 cells, and
999 permutations. The outcomes calibrate the identity-testing stage only: program
matches are stipulated and no SAE is trained. Five paired, common perturbation
offset checks assess rank invariance; they are not additional independent trials.

Example (run from the repository root):
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src \
    python src/calibrate_program_modulation_specificity.py \
    --out-dir data_program_modulation_audit/calibration

The output directory must not already exist. Partial results are preserved after
failure and must be reviewed; the script never overwrites or silently resumes.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
from scipy.stats import binomtest

import program_modulation_specificity as identity


GENERATOR_STREAM_SALT = 20261018


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _timestamp():
    return datetime.now(timezone.utc).isoformat()


def _write_json_new(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def _simulate(seed, parameters):
    n_pert_cells = parameters["n_perturbations"] * parameters["cells_per_perturbation"]
    n_cells = parameters["n_controls"] + n_pert_cells
    rng = np.random.default_rng(np.random.SeedSequence([seed, GENERATOR_STREAM_SALT, 1]))
    activations = rng.normal(size=(n_cells, parameters["n_programs"]))
    labels = np.array(["control"] * parameters["n_controls"] +
                      [f"P{p}" for p in range(parameters["n_perturbations"])
                       for _ in range(parameters["cells_per_perturbation"])])
    records = {f"program_{j}": {"feature": j, "enrich_p_bonf": 1e-6}
               for j in range(parameters["n_programs"])}
    return activations, labels, records


def _run_test(activations, labels, records, seed, parameters):
    return identity.program_modulation_specificity(
        activations, labels, records, seed=seed, n_shuffle=parameters["n_permutations"],
        min_cells=parameters["min_cells"], alpha=parameters["alpha"],
        match_alpha=parameters["match_alpha"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=Path("data_program_modulation_audit/calibration"))
    parser.add_argument("--n-simulations", type=int, default=200)
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--n-permutations", type=int, default=999)
    parser.add_argument("--progress-every", type=int, default=20)
    args = parser.parse_args(argv)
    if args.n_simulations < 1 or args.seed_start < 0 or args.n_permutations < 1 or args.progress_every < 1:
        parser.error("simulation count, permutation count, and progress interval must be positive; seed start must be nonnegative")
    if args.out_dir.exists():
        parser.error(f"refusing to overwrite existing output directory: {args.out_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=False)
    parameters = {
        "n_simulations": args.n_simulations, "seed_start": args.seed_start,
        "seed_end_inclusive": args.seed_start + args.n_simulations - 1,
        "n_permutations": args.n_permutations, "n_programs": 9, "n_distinct_features": 9,
        "n_controls": 100, "n_perturbations": 10, "cells_per_perturbation": 40,
        "min_cells": 20, "alpha": 0.05, "match_alpha": 0.05,
        "stipulated_marker_match_p_bonf": 1e-6, "generator": "independent standard normal float64 features for every cell",
        "generator_stream_salt": GENERATOR_STREAM_SALT,
        "generator_seed_sequence": "[simulation_seed, generator_stream_salt, 1]",
        "permutation_seed": "simulation_seed (separate generator stream)",
        "paired_generic_checks": min(5, args.n_simulations), "paired_generic_offset": 3.0,
        "strata": None, "multiplicity": "Holm over all nine planned programs",
        "outcome": "any Holm-adjusted conditional_heterogeneity rejection in one simulation",
        "confidence_interval": "two-sided exact Clopper-Pearson 95% interval for simulation familywise rejection frequency",
        "anti_conservatism_diagnostic": "one-sided exact binomial test against nominal familywise rate 0.05; no tuning or reruns",
    }
    script_path = Path(__file__).resolve()
    module_path = Path(identity.__file__).resolve()
    provenance = {
        "started_utc": _timestamp(), "argv": sys.argv if argv is None else [str(script_path), *argv],
        "working_directory": str(Path.cwd()), "parameters": parameters,
        "sources": {
            "calibration_generator": {"path": str(script_path), "sha256": _sha256(script_path)},
            "executed_identity_module": {"path": str(module_path), "sha256": _sha256(module_path)},
        },
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "scope": "Synthetic calibration of conditional cell-level identity testing with stipulated matches. Not biological evidence, marker-enrichment calibration, SAE-learning validation, or discovery.",
        "paired_generic_check_scope": "Paired rank-invariance checks share datasets with the null simulations and are not additional independent trials.",
    }
    _write_json_new(args.out_dir / "provenance.json", provenance)
    started = time.perf_counter()
    rejections = 0
    invariance_checks = 0
    with (args.out_dir / "simulation_results.jsonl").open("x", buffering=1) as results_file:
        for index, seed in enumerate(range(args.seed_start, args.seed_start + args.n_simulations)):
            sim_started = time.perf_counter()
            activations, labels, records = _simulate(seed, parameters)
            result = _run_test(activations, labels, records, seed, parameters)
            rejected = [name for name, rec in result["programs"].items() if rec["conditional_heterogeneity"]]
            row = {
                "simulation_seed": seed, "any_holm_rejection": bool(rejected), "rejected_programs": rejected,
                "n_permutations_executed": result["n_permutations"],
                "programs": {name: {key: rec[key] for key in (
                    "feature", "status", "empirical_p", "holm_p", "conditional_heterogeneity", "rank_eta_squared")}
                    for name, rec in result["programs"].items()},
            }
            if index < parameters["paired_generic_checks"]:
                shifted = activations.copy()
                shifted[labels != "control"] += parameters["paired_generic_offset"]
                shifted_result = _run_test(shifted, labels, records, seed, parameters)
                checks = [all(rec[key] == shifted_result["programs"][name][key]
                              for key in ("empirical_p", "holm_p", "conditional_heterogeneity", "rank_eta_squared"))
                          for name, rec in result["programs"].items()]
                row["paired_generic_offset_rank_invariant"] = all(checks)
                if not all(checks):
                    raise RuntimeError(f"generic-offset rank invariance failed at simulation seed {seed}; partial evidence preserved")
                invariance_checks += 1
            row["elapsed_seconds"] = round(time.perf_counter() - sim_started, 6)
            results_file.write(json.dumps(row, allow_nan=False) + "\n")
            rejections += bool(rejected)
            if (index + 1) % args.progress_every == 0 or index + 1 == args.n_simulations:
                print(json.dumps({"completed": index + 1, "planned": args.n_simulations,
                                  "familywise_rejections": rejections,
                                  "elapsed_seconds": round(time.perf_counter() - started, 3)}), flush=True)
    for source in provenance["sources"].values():
        if _sha256(Path(source["path"])) != source["sha256"]:
            raise RuntimeError(f"source changed during calibration: {source['path']}; partial evidence preserved")
    fit = binomtest(rejections, args.n_simulations, p=parameters["alpha"], alternative="greater")
    interval = binomtest(rejections, args.n_simulations).proportion_ci(confidence_level=0.95, method="exact")
    summary = {
        "completed_utc": _timestamp(), "n_simulations": args.n_simulations,
        "familywise_rejections": rejections, "familywise_rejection_frequency": rejections / args.n_simulations,
        "exact_95_percent_binomial_interval": [float(interval.low), float(interval.high)],
        "nominal_familywise_alpha": parameters["alpha"],
        "one_sided_binomial_p_for_rate_above_nominal": float(fit.pvalue),
        "anti_conservatism_signal_at_0_05": bool(fit.pvalue < 0.05),
        "paired_generic_offset_rank_invariance_checks_passed": invariance_checks,
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "sources": provenance["sources"], "parameters": parameters,
        "simulation_results_sha256": _sha256(args.out_dir / "simulation_results.jsonl"),
        "interpretation": "Finite synthetic null calibration only; neither a proof of general validity nor evidence for a biological effect or discovery.",
    }
    _write_json_new(args.out_dir / "summary.json", summary)
    print(json.dumps({"summary": str(args.out_dir / "summary.json"),
                      "familywise_rejections": rejections, "n_simulations": args.n_simulations,
                      "exact_95_percent_binomial_interval": summary["exact_95_percent_binomial_interval"],
                      "one_sided_binomial_p_for_rate_above_nominal": fit.pvalue}), flush=True)


if __name__ == "__main__":
    main()
