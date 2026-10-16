"""Glaser block calculation (steady vapour diffusion) for the walls rig's three model walls.

Scene: DIN 4108-3:2001, Annex A (simplified winter diffusion calculation), as transcribed by umwelt-online (a
third-party HTML transcription):
  winter (dew) period 1440 h; inside 20 C / 50 %, outside -10 C / 80 %; R_si 0.13, R_se 0.04 m2K/W;
  diffusion resistance Z = 1.5e6 * s_d  (m2 h Pa / kg), i.e. delta_0 = 1 / 1.5e6 kg / (m h Pa);
  saturation pressure (A.4.2): 288.68 (1.098 + t/100)^8.02 for t >= 0 C, 4.689 (1.486 + t/100)^12.30 below 0 C.
A later block climate (DIN 4108-3:2014; the same values for walls in 2018 and 2024: 90 days, inside 20 C / 50 %,
outside -5 C / 80 %, R_si 0.25, delta_0 = 2e-10 kg/(m s Pa), Magnus saturation pressure) is run only as a
boundary-condition sensitivity check. It is not a statement about which edition applies to which wall.

Walls: the rig's three model walls, read unchanged from ../../reference/result.json (same brick, same insulation,
same U-value; insulation outside, in the middle, or inside the brick). Vapour-resistance factors mu are not in the
parent model, so they are swept: brick 5, 10, 16, 50, 100 for each brick leaf independently (DIN 4108-4 solid and
perforated clay brick 5/10, ISO 10456 fired clay 16/10, DIN 4108-4 clinker 50/100); mineral wool 1.

Method: every layer is split into sub-layers of R <= 0.25 m2K/W (refine multiplies the count). Temperatures follow
the thermal resistances; saturation pressure is evaluated at every interface. The vapour-pressure line is the lower
convex hull of the points (s_d, p_sat) between the outdoor (0, p_e) and the indoor (s_d,total, p_i) air values -- the
Glaser construction: where the straight line would exceed saturation, the line is pulled down onto the saturation
curve and the touching interfaces are condensation planes. The condensation rate at a plane is the difference of the
diffusion flows into and out of it; times the dew period it gives the water-equivalent condensate (kg/m2).
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parents[1] / "reference" / "result.json"       # the rig's first study: the walls themselves

MU_BRICK = (5, 10, 16, 50, 100)
MU_BRICK_ORDINARY = (5, 10, 16)
MU_WOOL = 1.0
MU_NOMINAL = 10
REFINE = 16          # sub-layers per R <= 0.25 m2K/W slice; refine 1 is kept as the refinement check


def psat_din2001(t):
    """DIN 4108-3:2001 A.4.2."""
    return 288.68 * (1.098 + t / 100) ** 8.02 if t >= 0 else 4.689 * (1.486 + t / 100) ** 12.30


def psat_magnus(t):
    """Magnus form (water above 0 C, ice below), used by the later block climate."""
    if t >= 0:
        return 610.5 * math.exp(17.269 * t / (237.3 + t))
    return 610.5 * math.exp(21.875 * t / (265.5 + t))


BLOCKS = {
    "2001": dict(name="DIN 4108-3:2001 Annex A", ti=20.0, rhi=0.5, te=-10.0, rhe=0.8, hours=1440.0, rsi=0.13,
                 rse=0.04, d0=1 / 1.5e6, psat=psat_din2001),                    # d0: kg / (m h Pa)
    "2014": dict(name="DIN 4108-3:2014 block (sensitivity only)", ti=20.0, rhi=0.5, te=-5.0, rhe=0.8, hours=2160.0,
                 rsi=0.25, rse=0.04, d0=2e-10 * 3600, psat=psat_magnus),
}


def parent_model():
    return json.loads(PARENT.read_text(encoding="utf-8"))


def walls(mu_outer_brick=MU_NOMINAL, mu_inner_brick=MU_NOMINAL, mu_wool=MU_WOOL):
    """The parent's three walls, outside -> room, as (material, thickness m, lambda, mu). In the middle wall the two
    brick leaves take their own mu; in the single-leaf walls the brick takes the leaf it replaces (inside-insulated:
    the brick is the outer leaf; outside-insulated: the inner leaf)."""
    m = parent_model()
    lam = {k: v["lam"] for k, v in m["model"]["materials"].items()}
    out = {}
    for name, case in m["cases"].items():
        bricks = [i for i, (mat, _) in enumerate(case["layers"]) if mat == "brick"]
        lay = []
        for i, (mat, d) in enumerate(case["layers"]):
            if mat == "insulation":
                mu = mu_wool
            elif len(bricks) == 2:
                mu = mu_outer_brick if i == bricks[0] else mu_inner_brick
            else:
                mu = mu_outer_brick if name == "inside" else mu_inner_brick
            lay.append((mat, float(d), float(lam[mat]), float(mu)))
        out[name] = lay
    return out


def nodes(layers, rsi, rse, refine=1):
    """Interfaces from the outdoor surface to the indoor surface: x (mm), R from outdoor air, s_d, layer index."""
    x, R, sd, idx = [0.0], [rse], [0.0], [-1]
    for li, (_, d, lam, mu) in enumerate(layers):
        n = max(1, math.ceil((d / lam) / 0.25)) * refine
        for k in range(n):
            x.append(x[-1] + 1000 * d / n)
            R.append(R[-1] + d / lam / n)
            sd.append(sd[-1] + d * mu / n)
            idx.append(li)
    return x, R, sd, idx, R[-1] + rsi


def lower_hull(pts):
    hull = [pts[0]]
    for p in pts[1:]:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (y2 - y1) * (p[0] - x1) >= (p[1] - y1) * (x2 - x1):
                hull.pop()
            else:
                break
        hull.append(p)
    return hull


def envelope_by_slopes(pts):
    """Independent construction (for checking): from the current contact point, walk to the point that minimises the
    forward slope (ties: the farthest), until the indoor end is reached."""
    i, out = 0, [pts[0]]
    while i < len(pts) - 1:
        best, bj = None, None
        for j in range(i + 1, len(pts)):
            s = (pts[j][1] - pts[i][1]) / (pts[j][0] - pts[i][0])
            if best is None or s < best - 1e-15 or abs(s - best) <= 1e-15:
                best, bj = s, j
        out.append(pts[bj])
        i = bj
    return out


def condensate(hull, d0, hours):
    planes, rate = [], 0.0
    for k in range(1, len(hull) - 1):
        (xa, ya), (xb, yb), (xc, yc) = hull[k - 1], hull[k], hull[k + 1]
        g_in = d0 * (yc - yb) / (xc - xb)          # from the room side (kg / m2 h)
        g_out = d0 * (yb - ya) / (xb - xa)         # on to the outdoor side
        planes.append(dict(sd=xb, p=yb, g_in=g_in, g_out=g_out, g=g_in - g_out))
        rate += g_in - g_out
    return planes, rate * hours


def glaser(layers, block="2001", refine=1, ti=None, rhi=None, envelope=lower_hull):
    B = BLOCKS[block]
    ti = B["ti"] if ti is None else ti
    rhi = B["rhi"] if rhi is None else rhi
    x, R, sd, idx, Rt = nodes(layers, B["rsi"], B["rse"], refine)
    T = [B["te"] + (ti - B["te"]) * r / Rt for r in R]
    ps = [B["psat"](t) for t in T]
    pe, pi = B["rhe"] * B["psat"](B["te"]), rhi * B["psat"](ti)
    pts = [(sd[0], pe)] + [(sd[k], ps[k]) for k in range(1, len(sd) - 1)] + [(sd[-1], pi)]
    hull = envelope(pts)
    planes, M = condensate(hull, B["d0"], B["hours"])
    # vapour pressure at every interface: the hull, interpolated in s_d
    p = []
    for s in sd:
        for (xa, ya), (xb, yb) in zip(hull, hull[1:]):
            if xa - 1e-12 <= s <= xb + 1e-12:
                p.append(ya + (yb - ya) * (s - xa) / (xb - xa) if xb > xa else ya)
                break
    straight = [pe + (pi - pe) * s / sd[-1] for s in sd]
    clear = None if planes else min(ps[k] - straight[k] for k in range(len(sd)))
    for pl in planes:
        k = min(range(len(sd)), key=lambda j: abs(sd[j] - pl["sd"]))
        pl["x_mm"], pl["T"] = x[k], T[k]
    flux_room = B["d0"] * (hull[-1][1] - hull[-2][1]) / (hull[-1][0] - hull[-2][0])
    flux_out = B["d0"] * (hull[1][1] - hull[0][1]) / (hull[1][0] - hull[0][0])
    return dict(M=M, planes=planes, clear=clear, x_mm=x, T=T, ps=ps, p=p, straight=straight, sd=sd, layer=idx,
                pe=pe, pi=pi, flux_in=flux_room, flux_out=flux_out, T_cold_plane=min((pl["T"] for pl in planes),
                                                                                     default=None))


def benchmark():
    """Xella/Hebel technical handbook, Example 3 (https://technik.xella.de/media/books/handbuch/151/ and /152/):
    250 mm AAC (mu 5, lambda 0.14) + 9 mm bitumen (mu 20000, lambda 0.17), R_si 0.25 / R_se 0.04, -5/80 and 20/50,
    90 days: printed condensate 0.905 kg/m2 and evaporation capacity 1.002 kg/m2 (roof drying block: p_i = p_e =
    1200 Pa, wet zone 2000 Pa, 90 days)."""
    lay = [("bitumen", .009, .17, 20000.0), ("aac", .25, .14, 5.0)]
    r = glaser(lay, "2014")
    d0 = BLOCKS["2014"]["d0"]
    ev = d0 * ((2000.0 - 1200.0) / 1.25 + (2000.0 - 1200.0) / 180.0) * 2160.0
    return dict(M_c=r["M"], printed_M_c=0.905, M_ev=ev, printed_M_ev=1.002)


def sweep(block="2001", refine=1, mus=MU_BRICK):
    rows = []
    for mo, mi in itertools.product(mus, mus):
        W = walls(mo, mi)
        r = {k: glaser(v, block, refine) for k, v in W.items()}
        rows.append(dict(mu_outer=mo, mu_inner=mi, **{k: r[k]["M"] for k in r},
                         outside_clear=r["outside"]["clear"]))
    return rows


def summary(rows):
    inside_min = min(x["inside"] for x in rows)
    middle_max = max(x["middle"] for x in rows)
    return dict(cells=len(rows), inside=[inside_min, max(x["inside"] for x in rows)],
                middle=[min(x["middle"] for x in rows), middle_max], outside_max=max(x["outside"] for x in rows),
                outside_min_clearance_Pa=min(x["outside_clear"] for x in rows if x["outside_clear"] is not None),
                inside_over_middle_joint=inside_min / middle_max,
                reversal=any(x["middle"] >= x["inside"] or x["outside"] >= x["inside"] for x in rows))


def study():
    out = {"scene": BLOCKS["2001"]["name"], "benchmark": benchmark()}
    W = walls()
    U = {}
    for k, lay in W.items():
        Rt = sum(d / lam for _, d, lam, _ in lay) + BLOCKS["2001"]["rsi"] + BLOCKS["2001"]["rse"]
        U[k] = 1 / Rt
    out["U_W_m2K"] = U
    nominal = {}
    for k, lay in W.items():
        r = glaser(lay, "2001", REFINE)
        nominal[k] = dict(M=r["M"], planes=[{q: pl[q] for q in ("x_mm", "p", "T", "g_in", "g_out", "g")}
                                            for pl in r["planes"]],
                          flux_in=r["flux_in"], flux_out=r["flux_out"], clear=r["clear"],
                          profile={q: r[q] for q in ("x_mm", "T", "ps", "p", "straight", "layer")},
                          layers=[[m, d, lam, mu] for m, d, lam, mu in lay], pe=r["pe"], pi=r["pi"])
    out["nominal_mu10"] = nominal
    cells = {}
    for ref in (1, 16):
        rows = sweep("2001", ref)
        cells[str(ref)] = dict(rows=rows, all=summary(rows),
                               ordinary=summary([x for x in rows if x["mu_outer"] in MU_BRICK_ORDINARY
                                                 and x["mu_inner"] in MU_BRICK_ORDINARY]))
    out["sweep_2001"] = cells
    # independent envelope (minimum forward slopes) over every 2001 cell at both refinements
    dev = 0.0
    for ref in (1, 16):
        for mo, mi in itertools.product(MU_BRICK, MU_BRICK):
            for k, lay in walls(mo, mi).items():
                a = glaser(lay, "2001", ref)["M"]
                b = glaser(lay, "2001", ref, envelope=envelope_by_slopes)["M"]
                dev = max(dev, abs(a - b))
    out["independent_envelope_max_abs_diff"] = dev
    # numerical uncertainty within the prescribed calculation
    refine_dev = max(abs(a[k] - b[k]) / a[k] for a, b in zip(cells["1"]["rows"], cells["16"]["rows"])
                     for k in ("inside", "middle") if a[k] > 0)
    bm = out["benchmark"]
    bench_dev = abs(bm["M_c"] - bm["printed_M_c"]) / bm["printed_M_c"]
    magnus = {}
    for k, lay in W.items():
        B = dict(BLOCKS["2001"])
        BLOCKS["2001x"] = dict(B, psat=psat_magnus)
        magnus[k] = glaser(lay, "2001x", REFINE)["M"]
        del BLOCKS["2001x"]
    psat_dev = max(abs(magnus[k] - nominal[k]["M"]) / nominal[k]["M"] for k in ("inside", "middle"))
    u = refine_dev + bench_dev + psat_dev
    out["u"] = dict(refinement=refine_dev, benchmark=bench_dev, psat_formula=psat_dev, total=u,
                    need=max(0.10, 3 * u))
    # the registered claim (fallback accepted at the existence gate): inside exceeds both, at the joint adverse
    # corner of every brick mu cell (refined calculation); middle vs outside unresolved. The coarse sub-layering's
    # values are kept as the refinement check (u).
    s_all = [cells[str(REFINE)]["all"]]
    s_ord = [cells[str(REFINE)]["ordinary"]]
    inside_min = min(s["inside"][0] for s in s_all)
    middle_max = max(s["middle"][1] for s in s_all)
    ratio = inside_min / middle_max
    out["claims"] = dict(
        inside_over_either_min=ratio, inside_min=inside_min, middle_max=middle_max,
        gap=ratio - 1, clears=(ratio - 1) > out["u"]["need"],
        middle_range_all=[min(s["middle"][0] for s in s_all), middle_max],
        middle_range_ordinary=[min(s["middle"][0] for s in s_ord), max(s["middle"][1] for s in s_ord)],
        outside_max=max(s["outside_max"] for s in s_all),
        outside_min_clearance_Pa=min(s["outside_min_clearance_Pa"] for s in s_all),
        middle_vs_outside="unresolved" if min(s["middle"][0] for s in s_all) <= 0 else "middle above outside",
        reversal=any(s["reversal"] for s in s_all))
    # later block climate: a boundary-condition sensitivity check only
    rows14 = sweep("2014", REFINE)
    nom14 = {k: glaser(lay, "2014", REFINE)["M"] for k, lay in W.items()}
    out["sensitivity_2014_block"] = dict(summary=summary(rows14), nominal=nom14,
                                         same_order=all(x["inside"] > max(x["middle"], x["outside"]) for x in rows14))
    # controls
    uni = {k: glaser([(m, d, lam, 10.0) for m, d, lam, _ in lay], "2001", REFINE)["M"] for k, lay in W.items()}
    skin = glaser([("bitumen", .009, .17, 20000.0)] + W["outside"], "2001", REFINE)
    humid = {k: glaser(lay, "2001", REFINE, rhi=0.8)["M"] for k, lay in W.items()}
    out["controls"] = dict(uniform_mu10=uni, outside_wall_with_vapour_tight_outer_skin=skin["M"],
                           humid_interior_80pct=humid,
                           identical_inputs_identical_results=glaser(W["middle"], "2001", 1)["M"]
                           == glaser(W["middle"], "2001", 1)["M"])
    return out


if __name__ == "__main__":
    o = study()
    c = o["claims"]
    print("benchmark", o["benchmark"])
    print("nominal", {k: round(v["M"], 4) for k, v in o["nominal_mu10"].items()})
    print("claims", c)
    print("u", o["u"])
    print("2014 block", o["sensitivity_2014_block"]["summary"], o["sensitivity_2014_block"]["nominal"])
    print("controls", o["controls"])
    print("independent envelope", o["independent_envelope_max_abs_diff"])
