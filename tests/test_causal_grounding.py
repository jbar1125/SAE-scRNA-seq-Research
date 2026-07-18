"""Synthetic validation of the causal-grounding metric.

A world with KNOWN causal structure: 6 programs, each a regulator whose knockdown
suppresses exactly that program. The metric must (a) ground the 6 true regulators,
(b) reject a no-effect perturbation, (c) reject a non-specific GLOBAL perturbation that
suppresses everything (defeated by the match-confidence gate), and (d) collapse to
~chance under label shuffling. Run: python tests/test_causal_grounding.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import causal_grounding as cg

N_PROG, PROG_SIZE = 6, 10
N_GENES = N_PROG * PROG_SIZE
GENES = [f"g{i}" for i in range(N_GENES)]
RNG = np.random.default_rng(0)
NOISE = 0.08


def _cells(activity):
    """Given per-program activity (n_cells, N_PROG), build (activations, expression).
    activation f = program-f activity; gene g expression = its program's activity."""
    A = activity + RNG.normal(0, NOISE, activity.shape)
    X = activity[:, np.repeat(np.arange(N_PROG), PROG_SIZE)] + RNG.normal(0, NOISE, (activity.shape[0], N_GENES))
    return A, X


def _world():
    acts, exprs, labels = [], [], []

    def base(n):
        return RNG.uniform(0.8, 1.2, size=(n, N_PROG))

    A, X = _cells(base(400)); acts.append(A); exprs.append(X); labels += ["control"] * 400
    for f in range(N_PROG):                                  # 6 true regulators
        act = base(60); act[:, f] = RNG.uniform(0.0, 0.1, size=60)   # strong, specific
        A, X = _cells(act); acts.append(A); exprs.append(X)
        labels += [f"g{f * PROG_SIZE}"] * 60
    for j in range(40):                                      # background: no effect
        A, X = _cells(base(40)); acts.append(A); exprs.append(X)
        labels += [f"bg{j}"] * 40
    act = base(60) * 0.6; A, X = _cells(act)                 # global: all down moderately
    acts.append(A); exprs.append(X); labels += ["g59"] * 60  # non-specific
    A, X = _cells(base(60)); acts.append(A); exprs.append(X) # null: no change
    labels += ["g58"] * 60
    return np.vstack(acts), np.vstack(exprs), np.array(labels)


def _stress_world():
    """A world modelling the REAL essential-gene confound: 3 true regulators (each the
    sole suppressor of its own program) PLUS 5 non-specific 'stressor' perturbations that
    ALL suppress one SHARED program (a 'cell-health' axis), the way many essential-gene
    knockdowns collapse a common program. Unlike the uniform 'global' case, each stressor's
    downstream signature is CONCENTRATED on the shared program, so it confidently MATCHES
    and is SUPPRESSED-vs-control -- it passes every gate EXCEPT column specificity."""
    acts, exprs, labels = [], [], []

    def base(n):
        return RNG.uniform(0.8, 1.2, size=(n, N_PROG))

    A, X = _cells(base(400)); acts.append(A); exprs.append(X); labels += ["control"] * 400
    for f in range(3):                                       # 3 sole-regulator programs
        act = base(60); act[:, f] = RNG.uniform(0.0, 0.1, size=60)
        A, X = _cells(act); acts.append(A); exprs.append(X); labels += [f"g{f * PROG_SIZE}"] * 60
    for j in range(5):                                      # 5 stressors share program 5
        act = base(60); act[:, 5] = RNG.uniform(0.0, 0.1, size=60)
        A, X = _cells(act); acts.append(A); exprs.append(X); labels += [f"s{j}"] * 60
    return np.vstack(acts), np.vstack(exprs), np.array(labels)


def test_column_specificity():
    """Column specificity must reject the shared non-specific stressors while keeping the
    true regulators -- and WITHOUT the gate (alpha=1.0) the stressors would falsely ground,
    which is exactly the essential-gene confound seen on real Replogle data."""
    A, X, labels = _stress_world()
    true_regs = [f"g{f * PROG_SIZE}" for f in range(3)]
    stressors = [f"s{j}" for j in range(5)]
    tested = true_regs + stressors

    # NB: column specificity is DEFAULT OFF (v5 rollback) because on real Replogle it
    # removed a genuine regulator (GATA1); it is enabled here explicitly to test that it
    # still does what it was designed to do on the shared-feature confound it targets.
    on = {r["pert"]: r for r in cg.causal_grounding(A, X, labels, GENES, tested,
                                                    colspec_alpha=0.10)["per_perturbation"]}
    off = {r["pert"]: r for r in cg.causal_grounding(A, X, labels, GENES, tested,
                                                     colspec_alpha=1.0)["per_perturbation"]}
    # without the gate, the stressors pass match+floor+FDR -> they falsely ground (the bug)
    n_false_off = sum(off[s]["grounded"] for s in stressors)
    assert n_false_off >= 3, f"gate-off should reproduce the confound (>=3 false), got {n_false_off}"
    # with the gate ON, stressors are rejected, true regulators survive
    assert all(not on[s]["grounded"] for s in stressors), \
        f"column specificity must reject shared stressors: {[(s, on[s]['grounded']) for s in stressors]}"
    assert sum(on[r]["grounded"] for r in true_regs) >= 2, \
        f"true regulators must survive column specificity: {[(r, on[r]['grounded']) for r in true_regs]}"
    print(f"column specificity: gate-off false-grounds {n_false_off}/5 stressors -> "
          f"gate-on rejects all 5, keeps {sum(on[r]['grounded'] for r in true_regs)}/3 true regs")


def test_effect_size_control():
    """M1 core: the effect-size control must flag suppression that EXCEEDS the effect-size
    trend and must NOT flag suppression that is merely expected for a big effect. The crux
    assertion is decisive: two perturbations with the SAME suppression but different total
    effect size get OPPOSITE calls -- the modest-effect one is a genuine outlier, the
    high-effect one is on-trend."""
    rng = np.random.default_rng(0)
    n_null = 80
    eff_null = rng.uniform(1.0, 14.0, n_null)
    supp_null = 0.03 * eff_null + rng.normal(0, 0.008, n_null)   # suppression EXPLAINED by effect
    # true signal: MODEST effect, suppression far above the trend
    eff_true = np.array([4.0, 4.3, 4.6, 4.9, 5.2])
    supp_true = 0.03 * eff_true + 0.25
    # crux pair: identical suppression 0.39, but one at low effect (outlier) vs high effect (on-trend)
    S = 0.39
    eff_pair = np.array([4.0, 13.0]); supp_pair = np.array([S, S])
    eff = np.concatenate([eff_null, eff_true, eff_pair])
    supp = np.concatenate([supp_null, supp_true, supp_pair])

    p = cg._effect_matched_pvalue(supp, eff, k=20)
    ti = range(n_null, n_null + 5)
    lo_effect, hi_effect = n_null + 5, n_null + 6           # the crux pair indices

    assert all(p[i] < 0.10 for i in ti), f"true residual outliers must be flagged: {[p[i] for i in ti]}"
    assert p[lo_effect] < 0.10, f"same suppression at LOW effect must be flagged, got {p[lo_effect]:.3f}"
    assert p[hi_effect] > 0.20, f"same suppression at HIGH effect must NOT be flagged, got {p[hi_effect]:.3f}"
    # most effect-explained nulls are not flagged (calibration)
    frac_null_flagged = np.mean(p[:n_null] < 0.10)
    assert frac_null_flagged < 0.15, f"effect-explained nulls should rarely flag: {frac_null_flagged:.2f}"
    print(f"effect-size control: crux pair (supp={S}) -> low-effect p={p[lo_effect]:.3f} FLAGGED, "
          f"high-effect p={p[hi_effect]:.3f} not; null false-flag {frac_null_flagged:.0%}")


def main():
    A, X, labels = _world()
    true_regs = [f"g{f * PROG_SIZE}" for f in range(N_PROG)]
    tested = true_regs + ["g59", "g58"] + [f"bg{j}" for j in range(5)]

    res = cg.causal_grounding(A, X, labels, GENES, tested)
    calls = {r["pert"]: r["grounded"] for r in res["per_perturbation"]}

    n_true = sum(calls[r] for r in true_regs)
    assert n_true == N_PROG, f"expected 6 true regulators grounded, got {n_true}: {calls}"
    assert calls["g58"] is False, "null perturbation must not be grounded"
    assert calls["g59"] is False, "global (non-specific) perturbation must not be grounded"
    assert all(calls[f"bg{j}"] is False for j in range(5)), "background must not be grounded"

    # (d) calibration: a world with NO real causal structure (every "perturbation" is
    # drawn from the control distribution) must ground ~0 -> false positives controlled.
    nacts, nlab = [], []
    A0, _ = _cells(RNG.uniform(0.8, 1.2, size=(400, N_PROG)))
    nacts.append(A0); nlab += ["control"] * 400
    nexpr = [_cells(RNG.uniform(0.8, 1.2, size=(400, N_PROG)))[1]]
    for lab in tested:
        act = RNG.uniform(0.8, 1.2, size=(60, N_PROG)); a, x = _cells(act)
        nacts.append(a); nexpr.append(x); nlab += [lab] * 60
    res_n = cg.causal_grounding(np.vstack(nacts), np.vstack(nexpr), np.array(nlab), GENES, tested)
    assert res_n["n_grounded"] <= 1, f"null-world grounding should be ~0, got {res_n['n_grounded']}"

    print(f"true regulators grounded: {n_true}/{N_PROG}; null={calls['g58']}; "
          f"global={calls['g59']}; background all False={all(not calls[f'bg{j}'] for j in range(5))}")
    print(f"grounding rate (of {res['n_tested']} tested): {res['causal_grounding_rate']:.3f}")
    print(f"null-world grounding (calibration): {res_n['n_grounded']}/{res_n['n_tested']} (expect ~0)")

    test_column_specificity()
    test_effect_size_control()
    print("ALL CAUSAL-GROUNDING TESTS PASSED")


if __name__ == "__main__":
    main()
