"""Reproduce the reheating-after-time-away comparison from frozen sources.

    python reproduce.py        (from this folder; needs the walls rig's reproduce.py one level up)

Re-runs the three rooms after a 12 h and a 28-day absence, the pairwise rule with its input-sensitivity index u,
the shared-corner search, the swap time, the house-sized-heater rule, the heater-power sweep and every
validation check, then compares each published number with reference/result.json. Exits nonzero, naming the
quantity, if anything disagrees or a source changed.
"""
from __future__ import annotations

import hashlib
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "solver"))


def sha(p):
    return hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()   # line-ending independent


def close(a, b, rel=1e-6):
    return abs(a - b) <= rel * max(1.0, abs(b))


def main():
    ref = json.loads((HERE / "reference" / "result.json").read_text(encoding="utf-8"))
    for rel, want in ref["source_sha256"].items():
        if sha((HERE / rel).resolve()) != want:
            raise SystemExit(f"source changed: {rel}")
    import reheat
    with tempfile.TemporaryDirectory() as tmp:
        o = reheat.run(tmp)
    bad = []
    for leg, rooms in ref["recovery_min"].items():
        for k, v in rooms.items():
            if not close(o["legs"][leg]["recovery_min"][k], v):
                bad.append(f"recovery {leg} {k}")
            if not close(o["legs"][leg]["u"][k], ref["u"][leg][k]):
                bad.append(f"u {leg} {k}")
    if not close(o["swap"]["D_star_h"], ref["swap_h"]):
        bad.append("swap time")
    for a, b in zip(o["swap"]["tested_corner_band_h"], ref["swap_band_h"]):
        if not close(a, b):
            bad.append("swap corner band")
    for t, r_ in zip(o["heater_sweep"]["thresholds_12h"], ref["heater_thresholds_12h"]):
        if not close(t["P_W"], r_["P_W"], 1e-4) or t["direction"] != r_["direction"]:
            bad.append("heater sweep thresholds")
    for leg, pairs in ref["pairs"].items():
        for k, want in pairs.items():
            got = o["legs"][leg]["pairs"][k]
            if not (close(got["gap"], want["gap"]) and close(got["need"], want["need"]) and got["clears"] == want["clears"]
                    and got["leader"] == want["leader"]):
                bad.append(f"pair {leg} {k}")
    for leg, want in ref["corners"].items():
        got = o["corners"][leg]
        if not (close(got["shared_worst_ratio"], want["shared_worst_ratio"]) and got["flips"] == want["flips"]):
            bad.append(f"corners {leg}")
    alt = o["alternatives"]["sized_heater_20C"]
    if not close(alt["P_W"], ref["sized_heater"]["P_W"]):
        bad.append("sized heater power")
    for k, want in ref["sized_heater"].items():
        if k == "P_W":
            continue
        for q in ("outside_min", "inside_min", "ratio"):
            if not close(alt[k][q], want[q]):
                bad.append(f"sized heater {k} {q}")
    for leg, rooms in ref["recovery_furnished_min"].items():
        for k, v in rooms.items():
            if not close(o["legs"][leg]["recovery_furnished_min"][k], v):
                bad.append(f"furnished recovery {leg} {k}")
        if not close(o["headline_numbers"][leg]["furnished_ratio"], ref["furnished_ratio"][leg]):
            bad.append(f"furnished ratio {leg}")
    if o["claims"] != ref["claims"]:
        bad.append("claims")
    if not o["validation_pass"]:
        bad.append("validation")
    if bad:
        raise SystemExit("MISMATCH: " + "; ".join(bad))
    for leg, name in (("d12h", "away 12 h"), ("d28d", "away 28 days")):
        r = o["legs"][leg]["recovery_min"]
        print(f"{name:13s} air back at 20 C after: " + ", ".join(f"{k} {r[k]:.2f} min" for k in ("inside", "middle", "outside")))
    for leg in ("d12h", "d28d"):
        p = o["legs"][leg]["pairs"]["outside_vs_inside"]
        print(f"{leg}: outside vs inside gap {p['gap']:.0%} vs max(10%, 3u) = {p['need']:.0%} -> "
              + ("clears" if p["clears"] else "within uncertainty"))
    lo, hi = o["swap"]["tested_corner_band_h"]
    for leg, name in (("d12h", "12 h"), ("d28d", "28 days")):
        print(f"{name}: furnished (180 kJ/K) ratio x{o['headline_numbers'][leg]['furnished_ratio']:.1f}")
    print(f"sized heater {alt['P_W']:.0f} W, 12 h: outside {alt['d12h_furn0']['outside_min']:.1f} vs inside "
          f"{alt['d12h_furn0']['inside_min']:.1f} min (a tie)")
    print(f"swap: about {lo / 24:.0f}-{hi / 24:.0f} days away (recomputed over the tested shared corners, 1.5 kW fixed: "
          f"{lo:.1f}-{hi:.1f} h; nominal {o['swap']['D_star_h']:.1f} h)")
    th = o["heater_sweep"]["thresholds_12h"]
    print(f"12 h pair vs heater power: inside leads below {th[0]['P_W']:.0f} W, within uncertainty to {th[1]['P_W']:.0f} W, "
          f"outside leads and clears from {th[1]['P_W']:.0f} W to 3 kW")
    print("all validation checks pass; recomputed recovery times (empty and furnished), u, pairs, corners, swap band, "
          "sized-heater power and times, and heater thresholds match the frozen reference")


if __name__ == "__main__":
    main()
