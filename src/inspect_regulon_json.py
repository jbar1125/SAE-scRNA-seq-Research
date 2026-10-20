#!/usr/bin/env python3
"""Diagnose a regulon_recovery.json: WHICH gate zeroed the result, and is the null informative?

A 0-recovered result is only a real negative if the test could have fired. This prints, per TF,
the two gates separately (regulon-enrichment of the selected feature vs causal activation of it),
so an underpowered/blind null is distinguishable from a genuine one -- the distinction CLAUDE.md
requires. It also reports the panel-truncation loss (TFs dropped for too few in-panel targets),
which is the usual cause of a tiny tested pool.

Reads only the JSON(s) already written by src/regulon_recovery.py; no data or GPU needed.

Usage: python3 src/inspect_regulon_json.py causal_out/regulon_recovery_seed*.json
"""
import argparse
import json


def show(path):
    res = json.load(open(path))
    per = res["per_tf"]
    tested = [r for r in per if not r.get("low_power")]
    low = [r for r in per if r.get("low_power")]

    print("=" * 78)
    print(f"{path}: recovered {res['n_recovered']}/{res['n_tf_tested']}")
    sn = res.get("shuffle_null")
    if sn:
        print(f"  shuffle null: mean {sn['null_mean']:.2f} max {sn['null_max']} "
              f"p={sn['empirical_p']:.3f}")
        if res["n_recovered"] == 0 and sn["null_max"] == 0:
            print("  *** UNINFORMATIVE: real AND null are both 0 -- the test never fired. "
                  "This is NOT evidence of absence. ***")
    print(f"  low-power (dropped, too few in-panel targets): {len(low)}"
          + (f" -> {[r['tf'] for r in low][:12]}" if low else ""))

    if not tested:
        print("  no TFs tested at all -- panel truncation killed every regulon.")
        return

    print(f"\n  {'TF':<10}{'targets':>8}{'cells':>7}{'enrichQ':>10}{'auc':>7}{'actQ':>10}"
          f"  {'gate that failed':<22} regulon hits")
    n_enrich_ok = n_floor_ok = n_act_ok = 0
    for r in sorted(tested, key=lambda r: r.get("enrich_q", 1.0)):
        eq = r.get("enrich_q", 1.0)
        aq = r.get("mw_q", 1.0)
        auc = r.get("auc", float("nan"))
        e_ok, f_ok, a_ok = eq < 0.05, r.get("passes_floor", False), aq < 0.05
        n_enrich_ok += e_ok
        n_floor_ok += f_ok
        n_act_ok += a_ok
        if r.get("recovers"):
            fail = "-- RECOVERS --"
        elif not e_ok:
            fail = "enrichment (selection)"
        elif not f_ok:
            fail = f"auc floor ({auc:.2f})"
        elif not a_ok:
            fail = "activation p"
        else:
            fail = "?"
        hits = ", ".join(r.get("regulon_hits", [])[:5])
        print(f"  {r['tf']:<10}{r['n_targets_in_panel']:>8}{r['n_oe_cells']:>7}"
              f"{eq:>10.2g}{auc:>7.2f}{aq:>10.2g}  {fail:<22} {hits}")

    n = len(tested)
    print(f"\n  gate pass counts (of {n} tested): enrichment {n_enrich_ok} | "
          f"auc floor {n_floor_ok} | activation {n_act_ok}")
    print("  READ: enrichment 0 -> the SAE features do not concentrate curated regulons "
          "(selection step is the bottleneck; panel truncation and/or feature granularity).")
    print("        enrichment >0 but activation 0 -> regulon features exist but OE does NOT "
          "activate them (that IS a real causal negative).")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("json", nargs="+")
    args = ap.parse_args()
    for p in args.json:
        show(p)


if __name__ == "__main__":
    main()
