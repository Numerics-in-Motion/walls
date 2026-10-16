"""Reproduce the winter-condensate comparison from frozen sources.

    python reproduce.py        (from this folder; reads the rig's walls from ../reference/result.json)

Re-runs the outside check (Xella/Hebel handbook Example 3), the three walls under DIN 4108-3:2001 Annex A, the brick
vapour-resistance sweep (5, 10, 16, 50, 100 per brick leaf) at two sub-layer refinements, an independent envelope
construction, the later block climate (sensitivity only) and the controls, then compares every number with
reference/result.json. Exits nonzero, naming the quantity, if anything disagrees or a source changed.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "solver"))

SOURCES = ("solver/glaser.py", "../reference/result.json")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes().replace(b"\r\n", b"\n")).hexdigest()   # line-ending independent


def flat(o, pre=""):
    out = {}
    if isinstance(o, dict):
        for k, v in o.items():
            out.update(flat(v, f"{pre}/{k}"))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out.update(flat(v, f"{pre}[{i}]"))
    else:
        out[pre] = o
    return out


def same(a, b, rel=1e-9):
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str):
        return a == b
    return abs(a - b) <= rel * max(1.0, abs(b))


def main():
    ref = json.loads((HERE / "reference" / "result.json").read_text(encoding="utf-8"))
    for rel, want in ref["source_sha256"].items():
        if sha(HERE / rel) != want:
            raise SystemExit(f"source changed: {rel}")
    import glaser
    res = json.loads(json.dumps(glaser.study()))
    got, want = flat(res), flat({k: v for k, v in ref.items() if k != "source_sha256"})
    bad = sorted(k for k in set(got) | set(want) if k not in got or k not in want or not same(got[k], want[k]))
    if bad:
        raise SystemExit("MISMATCH: " + "; ".join(bad[:20]))
    b, c, u = res["benchmark"], res["claims"], res["u"]
    checks = {
        "outside check: condensate prints as 0.905": round(b["M_c"], 3) == b["printed_M_c"],
        "outside check: evaporation prints as 1.002": round(b["M_ev"], 3) == b["printed_M_ev"],
        "same U-value for all three walls": max(res["U_W_m2K"].values()) - min(res["U_W_m2K"].values()) < 1e-12,
        "independent envelope agrees": res["independent_envelope_max_abs_diff"] < 1e-12,
        "inside exceeds both by more than max(10 %, 3u)": c["clears"] and c["gap"] > u["need"],
        "no reversal in any cell": not c["reversal"],
        "outside-insulated wall collects nothing": c["outside_max"] == 0.0,
        "zero is not a blind detector (vapour-tight skin condenses)":
            res["controls"]["outside_wall_with_vapour_tight_outer_skin"] > 0,
    }
    for k, v in checks.items():
        print(("PASS  " if v else "FAIL  ") + k)
    if not all(checks.values()):
        raise SystemExit("a check failed")
    n = res["nominal_mu10"]
    print(f"Xella/Hebel Example 3: condensate {b['M_c']:.4f} kg/m2 (printed {b['printed_M_c']}), "
          f"evaporation {b['M_ev']:.4f} (printed {b['printed_M_ev']})")
    print(f"DIN 4108-3:2001, brick mu 10: insulation inside {n['inside']['M']:.3f}, middle {n['middle']['M']:.3f}, "
          f"outside {n['outside']['M']:.3f} kg/m2 of water-equivalent condensate after 1440 h")
    print(f"every brick mu 5-100 per leaf: inside {c['inside_min']:.3f}+, middle {c['middle_range_all'][0]:.3f}-"
          f"{c['middle_range_all'][1]:.3f}, outside {c['outside_max']:.0f} -> inside >= x{c['inside_over_either_min']:.4f} "
          f"either other (need {u['need']:.0%}); middle vs outside {c['middle_vs_outside']}")
    s14 = res["sensitivity_2014_block"]
    print(f"later block climate (sensitivity only): same order in every cell: {s14['same_order']}")
    print("Every number matches the frozen reference.")


if __name__ == "__main__":
    main()
