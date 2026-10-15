"""Reheating after time away: three model rooms on the walls rig. Study module.

The three parent study model rooms (same brick, same insulation, same U-value; insulation inside / middle / outside)
are held at 20 C air with 0 C outdoors (steady state, heater at the steady loss). The heating goes off for an
absence D; on return the parent study's 1.5 kW heater goes on until the air node FIRST returns to 20 C.

modified; its frozen reference is reproduced first.

    python reproduce.py   (from the folder above)

Solver: the rig's own finite-volume network (cell-centred, dx 1.25 mm, one air node) solved exactly for
piecewise-constant heater power by modal decomposition. Checks: the rig's own reproduce.run(); the parent study's frozen t18 by
the modal solve; the registered crossings by the rig's Crank-Nicolson scheme and by Laplace superposition of the
rig's thermal-quadrupole step response; exact energy books; massless-wall control.

Admin condition (10/14): heater power is a freely chosen value, so it is SWEPT (0.4-3.0 kW) for the 12 h
outside-vs-inside pair before any headline; P_min = lowest power at which that pair clears max(10 %, 3u).
The registered 1.5 kW is not changed.
"""
from __future__ import annotations

import os
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "1")

import contextlib
import importlib.util
import io
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, lu_factor, lu_solve
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
RIG_DIR = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parents[1] / "out"


def _load_rig():
    spec = importlib.util.spec_from_file_location("walls_rig_reproduce", RIG_DIR / "reproduce.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


W = _load_rig()                       # frozen parent study rig: only read and called
REF = W.R
ROOMS = ("inside", "middle", "outside")
BUILDS = {k: [(m, t) for m, t in REF["cases"][k]["layers"]] for k in ROOMS}   # outside -> room
PAIR = ("outside", "inside")

T_HOME = 20.0                 # held air temperature and endpoint (registered)
P_HEAT = 1500.0               # the parent study's heater (registered; stated condition)
D_LEGS = {"d12h": 12.0, "d28d": 672.0}     # registered absences, hours
FURN = 180e3                  # the parent study's furnished variant, J/K on the air node

NOM = dict(b_lam=0.77, b_rho=1700.0, b_c=800.0, i_lam=0.035, i_rho=25.0, i_c=1400.0,
           rsi=0.13, rse=0.04, c_furn=0.0, dx=0.00125)
INP = dict(                            # registered one-at-a-time half-ranges (lo, hi)
    b_lam=(0.60, 0.95), b_rho=(1500.0, 1900.0), b_c=(750.0, 900.0),
    i_lam=(0.032, 0.040), i_rho=(15.0, 35.0), i_c=(1030.0, 1450.0),
    rsi=(0.10, 0.17), rse=(0.02, 0.06),
    c_furn=(0.0, 180e3),               # structural: furniture on the air node (the parent study's own variant)
    dx=(0.00125, 0.0025),              # structural: grid
)
KEYS = tuple(INP)

HEATER_BASIS = [
    {"claim": "1,000 W and 1,500 W are among the most common portable electric heater ratings (US utility guidance).",
     "url": "https://otpco.com/newsroom/posts/otter-tail-power-company-offers-tips-for-buying-and-using-portable-electric-space-heaters",
     "source": "Otter Tail Power Company newsroom, 2019-01-03",
     "quote": "The most common are 1,000-watt and 1,500-watt units."},
    {"claim": "A 1,500 W heater is the source's worked example; it draws about 12.5 A.",
     "url": "https://www.army.mil/article/69664/electrical_space_heaters",
     "source": "U.S. Army Combat Readiness/Safety Center (A. Brown), 2011-11-21",
     "quote": "If one 1,500-watt heater uses approximately 12.5 amps"},
]


# ----------------------------------------------------------------------------------------------- solver
class Room:
    """The rig's FV network (outside -> room cells + air node), per m^2 of wall; exact modal solution."""

    def __init__(self, build, mat=None, dx=None, c_room=None, rsi=None, rse=None):
        mat = mat or W.MAT
        dx = dx or W.M["dx_m"]
        c_room = W.C_AIR if c_room is None else c_room
        rsi = W.RSI if rsi is None else rsi
        rse = W.RSE if rse is None else rse
        self.area = W.AREA
        lam, cap, mats = [], [], []
        for m, t in build:
            n = round(t / dx)
            lam += [mat[m]["lam"]] * n
            cap += [mat[m]["rho"] * mat[m]["c"] * dx] * n
            mats += [m] * n
        lam, cap = np.array(lam), np.array(cap)
        n = len(lam)
        g = 1 / (dx / (2 * lam[:-1]) + dx / (2 * lam[1:]))
        g_out, g_in = 1 / (rse + dx / (2 * lam[0])), 1 / (rsi + dx / (2 * lam[-1]))
        cv = np.r_[cap, c_room / self.area]
        K = np.zeros((n + 1, n + 1))
        K[0, 0] += g_out
        for i in range(n - 1):
            K[i, i] += g[i]; K[i + 1, i + 1] += g[i]; K[i, i + 1] -= g[i]; K[i + 1, i] -= g[i]
        K[n - 1, n - 1] += g_in; K[n, n] += g_in; K[n - 1, n] -= g_in; K[n, n - 1] -= g_in
        self.K, self.cv, self.n, self.ncell, self.dx = K, cv, n + 1, n, dx
        self.mats = mats
        self.e = np.zeros(n + 1); self.e[-1] = 1.0 / self.area     # W -> W/m2 on the air node
        self.lam, self.V = eigh(K, np.diag(cv))                     # K v = lam C v, V^T C V = I
        self.Kinv_e = np.linalg.solve(K, self.e)
        self.air_per_W = self.Kinv_e[-1]
        self.g_out = g_out

    def steady(self, T_air=T_HOME):
        return self.Kinv_e * (T_air / self.air_per_W)

    def evolve(self, T0, P, t):
        Tinf = self.Kinv_e * P
        a = self.V.T @ (self.cv * (T0 - Tinf))
        t = np.atleast_1d(np.asarray(t, float))
        return Tinf[None, :] + (np.exp(-np.outer(t, self.lam)) * a[None, :]) @ self.V.T

    def air(self, T0, P, t):
        Tinf = self.Kinv_e[-1] * P
        a = self.V.T @ (self.cv * (T0 - self.Kinv_e * P))
        t = np.atleast_1d(t)
        return Tinf + (np.exp(-np.outer(t, self.lam)) * a[None, :]) @ self.V[-1, :]

    def time_to(self, T0, P, target, tmax=3e7):
        """First time the air node reaches `target` (registered bracketing, m01)."""
        f = lambda t: float(self.air(T0, P, t)[0]) - target
        if f(0.0) >= 0:
            return 0.0
        for tm in (tmax, 1e9):
            ts = np.r_[0.0, np.geomspace(1.0, tm, 400)]
            vals = self.air(T0, P, ts) - target
            i = int(np.argmax(vals >= 0))
            if vals[i] >= 0:
                return brentq(f, ts[i - 1], ts[i], xtol=1e-6, rtol=1e-12)
        return math.inf

    def energy_books(self, T0, P, t):
        """Exact integrals over one constant-power phase, per m^2: (heat in, stored change, lost outdoors)."""
        Tinf = self.Kinv_e * P
        a = self.V.T @ (self.cv * (T0 - Tinf))
        integ = Tinf * t + self.V @ (a * (-np.expm1(-self.lam * t)) / self.lam)
        T1 = self.evolve(T0, P, t)[0]
        return P / self.area * t, float(self.cv @ (T1 - T0)), self.g_out * float(integ[0])


_cache: dict = {}


def room(name, p):
    key = (name,) + tuple(p[k] for k in KEYS)
    if key not in _cache:
        if len(_cache) > 400:
            _cache.clear()
        mat = {"brick": {"lam": p["b_lam"], "rho": p["b_rho"], "c": p["b_c"]},
               "insulation": {"lam": p["i_lam"], "rho": p["i_rho"], "c": p["i_c"]}}
        r = Room(BUILDS[name], mat=mat, dx=p["dx"], c_room=W.C_AIR + p["c_furn"], rsi=p["rsi"], rse=p["rse"])
        r._TD = {}
        _cache[key] = r
    return _cache[key]


def state_after(r, D_h):
    if D_h not in r._TD:
        r._TD[D_h] = r.evolve(r.steady(), 0.0, D_h * 3600.0)[0]
    return r._TD[D_h]


def recov(name, p, D_h, P=P_HEAT, target=T_HOME):
    """Minutes from heater-on until the air node first returns to `target`."""
    r = room(name, p)
    return r.time_to(state_after(r, D_h), P, target) / 60.0


# ----------------------------------------------------------------------------------------------- decision rule
def u_table(D_h, P=P_HEAT, names=ROOMS):
    base = {n: recov(n, NOM, D_h, P) for n in names}
    u2 = {n: 0.0 for n in names}
    swings = {}
    for k, (lo, hi) in INP.items():
        rl = {n: recov(n, {**NOM, k: lo}, D_h, P) for n in names}
        rh = {n: recov(n, {**NOM, k: hi}, D_h, P) for n in names}
        swings[k] = {n: [rl[n], rh[n]] for n in names}
        for n in names:
            u2[n] += (0.5 * abs(rh[n] - rl[n]) / base[n]) ** 2
    return base, {n: math.sqrt(v) for n, v in u2.items()}, swings


def pair_rule(base, u, a, b):
    lead, trail = (a, b) if base[a] < base[b] else (b, a)
    gap = (base[trail] - base[lead]) / base[lead]
    need = max(0.10, 3 * max(u[a], u[b]))
    return dict(pair=[a, b], leader=lead, trailer=trail, lead_min=base[lead], trail_min=base[trail],
                ratio=base[trail] / base[lead], gap=gap, need=need, u={a: u[a], b: u[b]},
                clears=bool(gap > need), margin=gap / need)


def swap_time(p, lo=6.0, hi=400.0):
    return brentq(lambda D: recov("outside", p, D) - recov("inside", p, D), lo, hi, xtol=1e-3)


def corners():
    """All 2^10 lo/hi corners, shared materials in both rooms: recovery at both legs and the swap time."""
    rows = []
    for bits in itertools.product((0, 1), repeat=len(KEYS)):
        p = dict(NOM)
        for k, b in zip(KEYS, bits):
            p[k] = INP[k][b]
        row = {"bits": bits, "c_furn": p["c_furn"]}
        for lab, D in D_LEGS.items():
            row[lab] = {n: recov(n, p, D) for n in PAIR}
        row["swap_h"] = swap_time(p)
        rows.append(row)
        _cache.clear()
    return rows


def corner_summary(rows, lab, lead, trail):
    worst = min(rows, key=lambda c: c[lab][trail] / c[lab][lead])
    out = {"leader": lead, "trailer": trail,
           "shared_worst_ratio": worst[lab][trail] / worst[lab][lead],
           "shared_worst_values_min": worst[lab],
           "shared_worst_corner": {k: INP[k][b] for k, b in zip(KEYS, worst["bits"])},
           "flips": {}, "independent": {}}
    for f in (0.0, FURN):
        sub = [c for c in rows if c["c_furn"] == f]
        out["flips"][f"furniture_{f/1e3:.0f}kJK"] = f"{sum(1 for c in sub if c[lab][trail] <= c[lab][lead])}/{len(sub)}"
        lw, tb = max(c[lab][lead] for c in sub), min(c[lab][trail] for c in sub)
        out["independent"][f"equal_furniture_{f/1e3:.0f}kJK"] = {
            "leader_worst_min": lw, "trailer_best_min": tb, "ratio": tb / lw, "holds": bool(tb > lw)}
    lw, tb = max(c[lab][lead] for c in rows), min(c[lab][trail] for c in rows)
    out["independent"]["furniture_free_per_room"] = {"leader_worst_min": lw, "trailer_best_min": tb,
                                                     "ratio": tb / lw, "holds": bool(tb > lw)}
    out["flips_total"] = sum(1 for c in rows if c[lab][trail] <= c[lab][lead])
    return out


# ----------------------------------------------------------------------------------------------- heater sweep
def heater_point(P):
    rec = {"P_W": P}
    for lab, D in D_LEGS.items():
        base, u, _ = u_table(D, P, PAIR)
        rec[lab] = pair_rule(base, u, *PAIR)
    return rec


def heater_sweep():
    grid = [round(x, 2) for x in np.arange(400.0, 3000.0 + 1e-6, 20.0)]
    pts = [heater_point(P) for P in grid]
    f = lambda P: (lambda r: r["gap"] - r["need"])(heater_point(P)["d12h"])
    thresholds = []
    for a, b in zip(pts[:-1], pts[1:]):
        fa, fb = a["d12h"]["gap"] - a["d12h"]["need"], b["d12h"]["gap"] - b["d12h"]["need"]
        if fa * fb < 0:
            Pc = brentq(f, a["P_W"], b["P_W"], xtol=0.01)
            thresholds.append({"P_W": Pc, "direction": "clears above" if fb > 0 else "fails above"})
    clears = [p["P_W"] for p in pts if p["d12h"]["clears"]]
    # lowest power from which the 12 h pair clears continuously up to 3 kW (None if it fails at 3 kW)
    P_min = None
    if pts[-1]["d12h"]["clears"]:
        last_fail = max([p["P_W"] for p in pts if not p["d12h"]["clears"]], default=None)
        if last_fail is None:
            P_min = grid[0]                     # clears over the whole sweep (threshold below the range)
        else:
            P_min = min(t["P_W"] for t in thresholds if t["P_W"] > last_fail)
    leaders = {p["P_W"]: p["d12h"]["leader"] for p in pts}
    # the reversal needs the 12 h pair to clear with OUTSIDE leading (28 d is led by inside at every power)
    low = [p["P_W"] for p in pts if p["d12h"]["clears"] and p["d12h"]["leader"] == "inside"]
    low_band = None
    if low:
        top = [t["P_W"] for t in thresholds if t["direction"] == "fails above" and t["P_W"] > max(low)]
        low_band = [min(low), min(top) if top else max(low)]
    rev_ok = all(p["d12h"]["leader"] == "outside" for p in pts if p["P_W"] >= (P_min or math.inf))
    return {
        "note": ("The 12 h pair clears in TWO power bands: a low band (from the bottom of the sweep) where the "
                 "INSIDE room leads, i.e. the same order as 28 d and no reversal, and the band from "
                 "P_min_reversal_W upward where the OUTSIDE room leads, i.e. the reversal. Between them (incl. "
                 "the sized 657 W heater) the pair is a tie under max(10 %, 3u)."),
        "P_min_reversal_W": P_min if rev_ok else None,
        "low_band_clears_inside_leads_no_reversal_W": low_band,
        "steady_loss_W": None,
        "grid_W": grid,
        "points": [{"P_W": p["P_W"],
                    "d12h": {k: p["d12h"][k] for k in ("leader", "lead_min", "trail_min", "ratio", "gap", "need",
                                                       "u", "clears")},
                    "d28d": {k: p["d28d"][k] for k in ("leader", "lead_min", "trail_min", "ratio", "gap", "need",
                                                       "u", "clears")}} for p in pts],
        "thresholds_12h": thresholds,
        "P_min_W_12h_clears_continuously_to_3kW": P_min,
        "first_grid_point_clearing_12h_W": min(clears) if clears else None,
        "d28d_clears_all_grid_points": all(p["d28d"]["clears"] for p in pts),
        "d28d_min_margin": min(p["d28d"]["margin"] for p in pts),
        "d12h_leader_outside_all_grid_points": all(v == "outside" for v in leaders.values()),
        "d12h_powers_where_inside_leads_W": [P for P, v in leaders.items() if v != "outside"],
    }


# ----------------------------------------------------------------------------------------------- validation
def cn_crossing(r, D_h, P=P_HEAT):
    """The rig's scheme (4 implicit quarter steps, then CN), 60 s away / 1 s on return; with CN energy books."""
    K, cv = r.K, r.cv
    cache = {}
    books = {"in": 0.0, "lost": 0.0}

    def step(T, Pw, h, th):
        if (h, th) not in cache:
            cache[(h, th)] = lu_factor(np.diag(cv) + th * h * K)
        new = lu_solve(cache[(h, th)], cv * T - (1 - th) * h * (K @ T) + h * r.e * Pw)
        books["lost"] += r.g_out * (th * new[0] + (1 - th) * T[0]) * h
        books["in"] += r.e[-1] * Pw * h
        return new
    T = T0 = r.steady()
    for _ in range(4):
        T = step(T, 0.0, 15.0, 1.0)
    t = 60.0
    while t < D_h * 3600 - 1e-9:
        h = min(60.0, D_h * 3600 - t); T = step(T, 0.0, h, 0.5); t += h
    tt = 0.0
    seq = itertools.chain([(0.25, 1.0)] * 4, itertools.repeat((1.0, 0.5)))
    for h, th in seq:
        Tn = step(T, P, h, th)
        if T[-1] < T_HOME <= Tn[-1]:
            resid = books["in"] - float(cv @ (Tn - T0)) - books["lost"]
            return tt + h * (T_HOME - T[-1]) / (Tn[-1] - T[-1]), resid / books["lost"]
        T = Tn; tt += h


def laplace_crossing(name, r, D_h, guess):
    b = BUILDS[name]
    a = (T_HOME / r.air_per_W) / W.P_W
    g = lambda t: W.laplace_air(b, t) if t > 0 else 0.0
    f = lambda tau: -a * g(D_h * 3600 + tau) + g(tau)
    return brentq(f, 0.8 * guess, 1.25 * guess, xtol=1e-4)


def modal_energy_residual(r, D_h, tau_s):
    T0 = r.steady()
    hin1, st1, lo1 = r.energy_books(T0, 0.0, D_h * 3600)
    TD = r.evolve(T0, 0.0, D_h * 3600)[0]
    hin2, st2, lo2 = r.energy_books(TD, P_HEAT, tau_s)
    return (hin1 + hin2 - st1 - st2 - lo1 - lo2) / (lo1 + lo2)


def validation():
    v = {}
    # (0) the rig's own reproduce.run(): CN + Laplace + energy books vs its frozen reference (tol 1e-9)
    t0 = time.time()
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            W.run()
        ok = "EVERYTHING REPRODUCES" in buf.getvalue()
        msg = buf.getvalue().strip().splitlines()[-1]
    except AssertionError as exc:          # the rig names the failing quantity
        ok, msg = False, str(exc)
    v["rig_reproduce_run"] = {"value": msg, "seconds": time.time() - t0, "tol": REF["tolerance_relative"],
                              "pass": bool(ok)}
    # (1) parent study frozen t18 by the modal solve of the rig's FV system (cold start, 1.5 kW, 18 C), empty + furnished
    worst = 0.0; rows = {}
    for k in ROOMS:
        for furn, ref in ((0.0, REF["cases"][k]["t18_s"]), (FURN, REF["furnished_t18_s"][k])):
            r = Room(BUILDS[k], c_room=W.C_AIR + furn)
            t = r.time_to(np.zeros(r.n), W.P_W, W.TARGET)
            rel = t / ref - 1
            rows[f"{k}_{'furnished' if furn else 'empty'}"] = {"modal_s": t, "frozen_s": ref, "rel": rel}
            worst = max(worst, abs(rel))
    v["frozen_parent study_t18_modal"] = {"value": worst, "tol": 1e-6, "pass": bool(worst <= 1e-6), "cases": rows}
    # (2) registered crossings: modal vs rig CN vs Laplace superposition; (3) energy books
    cn_w = lap_w = en_w = cn_en_w = 0.0; rows = {}
    for lab, D in D_LEGS.items():
        for n in ROOMS:
            r = Room(BUILDS[n])
            m = r.time_to(state_after_plain(r, D), P_HEAT, T_HOME)
            c, cres = cn_crossing(r, D)
            l = laplace_crossing(n, r, D, m)
            e = modal_energy_residual(r, D, m)
            rows[f"{lab}_{n}"] = {"modal_min": m / 60, "cn_min": c / 60, "cn_rel": c / m - 1,
                                  "laplace_min": l / 60, "laplace_rel": l / m - 1,
                                  "modal_energy_resid": e, "cn_energy_resid": cres}
            cn_w, lap_w = max(cn_w, abs(c / m - 1)), max(lap_w, abs(l / m - 1))
            en_w, cn_en_w = max(en_w, abs(e)), max(cn_en_w, abs(cres))
    v["registered_cn_vs_modal"] = {"value": cn_w, "tol": 1e-4, "pass": bool(cn_w <= 1e-4), "cases": rows}
    v["registered_laplace_vs_modal"] = {"value": lap_w, "tol": 1e-3, "pass": bool(lap_w <= 1e-3)}
    v["energy_books_modal_exact"] = {"value": en_w, "tol": 1e-9, "pass": bool(en_w <= 1e-9),
                                     "note": "(heat in - stored change - lost outdoors) / lost, steady start to crossing"}
    v["energy_books_cn"] = {"value": cn_en_w, "tol": 1e-9, "pass": bool(cn_en_w <= 1e-9)}
    # (4) massless-wall control: identical recovery in all three rooms
    mat = {m: dict(W.MAT[m], rho=1e-6) for m in W.MAT}
    ctl = {}
    spread = 0.0
    for lab, D in D_LEGS.items():
        vals = {}
        for n in ROOMS:
            r = Room(BUILDS[n], mat=mat)
            vals[n] = r.time_to(state_after_plain(r, D), P_HEAT, T_HOME) / 60
        ctl[lab] = vals
        spread = max(spread, (max(vals.values()) - min(vals.values())) / min(vals.values()))
    v["massless_control_identical"] = {"value": spread, "tol": 1e-6, "pass": bool(spread <= 1e-6), "minutes": ctl}
    return v


def state_after_plain(r, D_h):
    return r.evolve(r.steady(), 0.0, D_h * 3600.0)[0]


# ----------------------------------------------------------------------------------------------- frames
def frames():
    arr = {}
    rooms = {n: Room(BUILDS[n]) for n in ROOMS}
    r0 = rooms["inside"]
    ncell = r0.ncell
    x_mm = (np.arange(ncell) + 0.5) * r0.dx * 1e3
    arr["rooms"] = np.array(ROOMS)
    arr["x_mm"] = x_mm
    arr["cell_material"] = np.array([[0 if m == "brick" else 1 for m in rooms[n].mats] for n in ROOMS], np.int8)
    arr["material_names"] = np.array(["brick", "insulation"])
    for n in ROOMS:
        edges = np.r_[0.0, np.cumsum([t for _, t in BUILDS[n]])] * 1e3
        arr[f"layer_bounds_mm_{n}"] = edges
        arr[f"layer_materials_{n}"] = np.array([m for m, _ in BUILDS[n]])
    arr["steady_air_C"] = np.full(3, T_HOME)
    arr["steady_wall_C"] = np.array([rooms[n].steady()[:-1] for n in ROOMS], np.float32)
    arr["steady_heater_W"] = np.array([T_HOME / rooms[n].air_per_W for n in ROOMS])
    arr["heater_W"] = np.array(P_HEAT)
    for lab, D in D_LEGS.items():
        Ds = D * 3600.0
        ta = np.r_[0.0, np.geomspace(60.0, Ds, 240)]
        ab = np.array([rooms[n].evolve(rooms[n].steady(), 0.0, ta) for n in ROOMS])   # (3, 241, ncell+1)
        arr[f"{lab}_D_s"] = np.array(Ds)
        arr[f"{lab}_absence_t_s"] = ta
        arr[f"{lab}_absence_air_C"] = ab[:, :, -1]
        arr[f"{lab}_absence_wall_C"] = ab[:, :, :-1].astype(np.float32)
        TD = {n: ab[i, -1] for i, n in enumerate(ROOMS)}
        tau = np.array([rooms[n].time_to(TD[n], P_HEAT, T_HOME) for n in ROOMS])
        arr[f"{lab}_crossing_s"] = tau
        own_t = np.array([np.linspace(0.0, tau[i], 200) for i in range(3)])
        own = np.array([rooms[n].evolve(TD[n], P_HEAT, own_t[i]) for i, n in enumerate(ROOMS)])
        arr[f"{lab}_reheat_own_t_s"] = own_t
        arr[f"{lab}_reheat_own_air_C"] = own[:, :, -1]
        arr[f"{lab}_reheat_own_wall_C"] = own[:, :, :-1].astype(np.float32)
        ts = np.unique(np.r_[own_t.ravel(), np.linspace(0.0, tau.max(), 200)])
        sh = np.array([rooms[n].evolve(TD[n], P_HEAT, np.minimum(ts, tau[i])) for i, n in enumerate(ROOMS)])
        arr[f"{lab}_reheat_shared_t_s"] = ts
        arr[f"{lab}_reheat_shared_air_C"] = sh[:, :, -1]
        arr[f"{lab}_reheat_shared_wall_C"] = sh[:, :, :-1].astype(np.float32)
        arr[f"{lab}_reheat_shared_done"] = ts[None, :] >= tau[:, None]
    return arr


# ----------------------------------------------------------------------------------------------- run
def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    return o


def run(out_dir=OUT_DIR):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    res = {"title": "Reheating after time away (walls rig)",
           "model": {"rig": "walls/reproduce.py (unchanged)", "materials": W.MAT,
                     "R_si": W.RSI, "R_se": W.RSE, "area_m2": W.AREA, "C_air_J_K": W.C_AIR, "dx_m": W.M["dx_m"],
                     "builds_outside_to_room": BUILDS, "T_home_C": T_HOME, "T_out_C": 0.0, "heater_W": P_HEAT,
                     "absences_h": D_LEGS, "endpoint": "air node first returns to 20 C after heater-on",
                     "U_W_m2K": REF["U_W_m2K"]},
           "inputs": {"nominal": NOM, "half_ranges": INP}}
    r0 = Room(BUILDS["inside"])
    P_ss = T_HOME / r0.air_per_W
    res["model"]["steady_loss_W"] = P_ss

    # registered legs: recovery per room, empty and furnished, u, pairwise table
    legs = {}
    for lab, D in D_LEGS.items():
        base, u, sw = u_table(D)
        furn = {n: recov(n, {**NOM, "c_furn": FURN}, D) for n in ROOMS}
        air_at_return = {n: float(state_after_plain(Room(BUILDS[n]), D)[-1]) for n in ROOMS}
        pairs = {f"{a}_vs_{b}": pair_rule(base, u, a, b)
                 for a, b in (("outside", "inside"), ("middle", "outside"), ("middle", "inside"))}
        for k, pr in pairs.items():
            a, b = pr["pair"]
            pr["furnished_min"] = {a: furn[a], b: furn[b]}
            pr["furnished_ratio"] = max(furn[a], furn[b]) / min(furn[a], furn[b])
            pr["furnished_leader"] = a if furn[a] < furn[b] else b
        legs[lab] = {"D_h": D, "recovery_min": base, "recovery_furnished_min": furn, "u": u,
                     "air_at_return_C": air_at_return, "swings_min": sw, "pairs": pairs}
    res["legs"] = legs

    # adversarial corners and swap time
    t0 = time.time()
    rows = corners()
    res["corners"] = {lab: corner_summary(rows, lab, legs[lab]["pairs"]["outside_vs_inside"]["leader"],
                                          legs[lab]["pairs"]["outside_vs_inside"]["trailer"]) for lab in D_LEGS}
    d0 = swap_time(NOM)
    uu = 0.0; oat = {}
    for k, (lo, hi) in INP.items():
        dl, dh = swap_time({**NOM, k: lo}), swap_time({**NOM, k: hi})
        uu += (0.5 * abs(dh - dl) / d0) ** 2; oat[k] = [dl, dh]
    sw_all = [c["swap_h"] for c in rows]
    res["swap"] = {"D_star_h": d0, "u": math.sqrt(uu), "one_at_a_time_h": oat,
                   "one_at_a_time_range_h": [min(min(v) for v in oat.values()), max(max(v) for v in oat.values())],
                   "tested_corner_band_h": [min(sw_all), max(sw_all)],
                   "tested_corner_band_d": [min(sw_all) / 24, max(sw_all) / 24],
                   "label": "tested shared corners (1024), 1.5 kW heater fixed; not a confidence interval",
                   "legs_outside_band": bool(12.0 < min(sw_all) and 672.0 > max(sw_all)),
                   "corner_seconds": time.time() - t0}

    # alternative rules and endpoint
    P2 = 2 * P_ss
    alt = {}
    for name, P, tgt in (("rule_parent study_1500W_20C", P_HEAT, T_HOME), ("sized_heater_20C", P2, T_HOME),
                         ("1500W_18C", P_HEAT, 18.0), ("sized_heater_18C", P2, 18.0)):
        alt[name] = {"P_W": P, "target_C": tgt}
        for lab, D in D_LEGS.items():
            for f in (0.0, FURN):
                o = recov("outside", {**NOM, "c_furn": f}, D, P, tgt)
                i = recov("inside", {**NOM, "c_furn": f}, D, P, tgt)
                alt[name][f"{lab}_furn{f/1e3:.0f}"] = {"outside_min": o, "inside_min": i,
                                                       "ratio": max(o, i) / min(o, i),
                                                       "leader": "outside" if o < i else "inside"}
    sized = heater_point(P2)
    alt["sized_heater_20C"]["pair_rule_12h"] = sized["d12h"]
    alt["sized_heater_20C"]["pair_rule_28d"] = sized["d28d"]
    res["alternatives"] = alt

    # admin condition 10/14: heater sweep
    t0 = time.time()
    hs = heater_sweep()
    hs["seconds"] = time.time() - t0
    hs["steady_loss_W"] = P_ss
    at15 = next(p for p in hs["points"] if abs(p["P_W"] - P_HEAT) < 1e-9)
    hs["at_1500W"] = at15
    res["heater_sweep"] = hs
    res["heater_basis"] = HEATER_BASIS

    # validation (kill line K2)
    t0 = time.time()
    val = validation()
    val["seconds"] = time.time() - t0
    res["validation"] = val
    res["validation_pass"] = all(v["pass"] for k, v in val.items() if isinstance(v, dict) and "pass" in v)

    # claims
    L12, L28 = legs["d12h"], legs["d28d"]
    oi12, oi28 = L12["pairs"]["outside_vs_inside"], L28["pairs"]["outside_vs_inside"]
    no_flips = all(res["corners"][lab]["flips_total"] == 0 for lab in D_LEGS)
    r28 = L28["recovery_min"]
    res["claims"] = {
        "pair_outside_inside_12h_clears": oi12["clears"],
        "pair_outside_inside_28d_clears": oi28["clears"],
        "reversal_claimed": bool(oi12["clears"] and oi28["clears"] and oi12["leader"] != oi28["leader"]
                                 and no_flips and res["validation_pass"]),
        "middle_12h_within_uncertainty": bool(not L12["pairs"]["middle_vs_outside"]["clears"]
                                              and not L12["pairs"]["middle_vs_inside"]["clears"]),
        "order_28d_inside_middle_outside": bool(r28["inside"] < r28["middle"] < r28["outside"]
                                                and all(p["clears"] for p in L28["pairs"].values())),
        "shared_corner_flips_zero": no_flips,
        "sized_heater_12h_is_tie_no_reversal_claimed_under_it": bool(not sized["d12h"]["clears"]),
        "heater_sweep_1500W_clears_12h": bool(at15["d12h"]["clears"]),
        "heater_sweep_28d_clears_everywhere": hs["d28d_clears_all_grid_points"],
        "heater_sweep_reversal_holds_from_P_min_to_3kW": hs["P_min_reversal_W"] is not None,
    }
    res["headline_numbers"] = {
        "d12h": {"outside_min": round(L12["recovery_min"]["outside"], 2),
                 "inside_min": round(L12["recovery_min"]["inside"], 1),
                 "ratio": round(oi12["ratio"]), "ratio_exact": oi12["ratio"],
                 "furnished_ratio": oi12["furnished_ratio"]},
        "d28d": {"inside_min": round(r28["inside"], 1), "outside_min": round(r28["outside"]),
                 "ratio": round(oi28["ratio"]), "ratio_exact": oi28["ratio"],
                 "furnished_ratio": oi28["furnished_ratio"]},
        "swap_band_days": [round(res["swap"]["tested_corner_band_d"][0]), round(res["swap"]["tested_corner_band_d"][1])],
        "swap_band_h": res["swap"]["tested_corner_band_h"],
        "heater_W": P_HEAT,
        "P_min_reversal_W": hs["P_min_reversal_W"],
        "heater_1500W_over_P_min": P_HEAT / hs["P_min_reversal_W"] if hs["P_min_reversal_W"] else None,
    }

    t0 = time.time()
    fr = frames()
    np.savez_compressed(out_dir / "frames.npz", **fr)
    res["frames"] = {"file": "frames.npz", "seconds": time.time() - t0,
                     "keys": {k: [list(v.shape), str(v.dtype)] for k, v in fr.items()}}
    res["run_seconds"] = time.time() - t_start
    (out_dir / "results.json").write_text(json.dumps(_clean(res), indent=1, ensure_ascii=False), encoding="utf-8")
    return res


if __name__ == "__main__":
    r = run(OUT_DIR)
    h = r["headline_numbers"]
    print("12 h: outside %.4f / inside %.4f min" % (r["legs"]["d12h"]["recovery_min"]["outside"],
                                                    r["legs"]["d12h"]["recovery_min"]["inside"]))
    print("28 d: inside %.4f / outside %.3f min" % (r["legs"]["d28d"]["recovery_min"]["inside"],
                                                    r["legs"]["d28d"]["recovery_min"]["outside"]))
    print("headline", json.dumps(_clean(h)))
    print("claims", json.dumps(r["claims"]))
    print("validation_pass", r["validation_pass"], {k: v.get("value") for k, v in r["validation"].items()
                                                     if isinstance(v, dict)})
    hs = r["heater_sweep"]
    print("P_min_reversal_W", hs["P_min_reversal_W"], "low band (inside leads)",
          hs["low_band_clears_inside_leads_no_reversal_W"], hs["thresholds_12h"])
    print("run %.1f s" % r["run_seconds"])
