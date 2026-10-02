"""Same brick, same insulation, same U-value: which model room's air reaches 18 C first?

    python reproduce.py

Solves the three walls with finite volumes, checks each finish time against an independent
Laplace-domain (thermal quadrupole) solution, closes the energy books, and compares every number
with reference/result.json. Exits nonzero, naming the quantity, if anything disagrees.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

R = json.loads((Path(__file__).parent / "reference" / "result.json").read_text())
M = R["model"]
MAT = M["materials"]
AREA, C_AIR, P_W = M["area_m2"], M["C_air_J_K"], M["heater_W"]
RSI, RSE, TARGET = M["R_si"], M["R_se"], M["target_C"]


def check(label, actual, expected, tol):
    if not math.isfinite(actual) or abs(actual - expected) > tol * max(1.0, abs(expected)):
        raise AssertionError("%s: %r != %r" % (label, actual, expected))


def r_total(build):
    return RSE + sum(t / MAT[m]["lam"] for m, t in build) + RSI


def finite_volumes(build, dx, dt, c_room=C_AIR):
    """Outside -> room cells, plus the room air node; Crank-Nicolson after 4 implicit substeps."""
    lam, cap = [], []
    for m, t in build:
        n = round(t / dx)
        lam += [MAT[m]["lam"]] * n
        cap += [MAT[m]["rho"] * MAT[m]["c"] * dx] * n
    lam, cap = np.array(lam), np.array(cap)
    n = len(lam)
    g = 1 / (dx / (2 * lam[:-1]) + dx / (2 * lam[1:]))
    g_out, g_in = 1 / (RSE + dx / (2 * lam[0])), 1 / (RSI + dx / (2 * lam[-1]))
    cv = np.r_[cap, c_room / AREA]
    d = np.zeros(n + 1)
    off = np.zeros(n)
    d[0] += g_out
    d[:n - 1] += g
    d[1:n] += g
    off[:n - 1] -= g
    d[n - 1] += g_in
    d[n] += g_in
    off[n - 1] -= g_in
    src = np.zeros(n + 1)
    src[n] = P_W / AREA

    def k(T):
        out = d * T
        out[:-1] += off * T[1:]
        out[1:] += off * T[:-1]
        return out
    T = np.zeros(n + 1)
    t, heat, lost, k_step = 0.0, 0.0, 0.0, 0
    while True:
        for h, th in ([(dt / 4, 1.0)] * 4 if k_step == 0 else [(dt, 0.5)]):
            ab = np.zeros((3, n + 1))
            ab[0, 1:], ab[1], ab[2, :-1] = th * h * off, cv + th * h * d, th * h * off
            new = solve_banded((1, 1), ab, cv * T - (1 - th) * h * k(T) + h * src)
            lost += g_out * (th * new[0] + (1 - th) * T[0]) * h
            heat += src[n] * h
            if T[-1] < TARGET <= new[-1]:
                tc = t + h * (TARGET - T[-1]) / (new[-1] - T[-1])
                stored = float(np.dot(cv, new))
                return tc, (heat - stored - lost) / heat
            T, t = new, t + h
        k_step += 1


def laplace_air(build, t):
    """Room-air rise from layer transfer matrices, inverted numerically (Talbot)."""
    def f(s):
        A = mp.matrix([[1, RSI], [0, 1]])
        for m, th in reversed(build):
            a = MAT[m]["lam"] / (MAT[m]["rho"] * MAT[m]["c"])
            q = mp.sqrt(s / a)
            A = A * mp.matrix([[mp.cosh(q * th), mp.sinh(q * th) / (MAT[m]["lam"] * q)],
                               [MAT[m]["lam"] * q * mp.sinh(q * th), mp.cosh(q * th)]])
        A = A * mp.matrix([[1, RSE], [0, 1]])
        return (mp.mpf(P_W) / AREA) / (s * (mp.mpf(C_AIR) / AREA * s + A[1, 1] / A[0, 1]))
    return float(mp.invertlaplace(f, t, method="talbot"))


def run():
    times = {}
    u0 = None
    for name, case in R["cases"].items():
        build = [(m, t) for m, t in case["layers"]]
        u = 1 / r_total(build)
        u0 = u if u0 is None else u0
        check("%s U-value" % name, u, u0, 1e-12)
        tc, resid = finite_volumes(build, M["dx_m"], M["dt_s"])
        check("%s energy closure" % name, resid, 0.0, 1e-9)
        lap = brentq(lambda t: laplace_air(build, t) - TARGET, 0.9 * tc, 1.1 * tc, xtol=1e-6)
        check("%s finite volumes vs Laplace" % name, tc / lap, 1.0, 1e-3)
        check("%s t18_s" % name, tc, case["t18_s"], R["tolerance_relative"])
        times[name] = tc
        print("%-8s %-46s U=%.4f W/m2K  air reaches 18 C in %s"
              % (name, " + ".join("%s %d mm" % (m, round(1e3 * t)) for m, t in build), u,
                 ("%.1f min" % (tc / 60)) if tc < 3600 else ("%.2f h" % (tc / 3600))))
    if not times["inside"] < times["middle"] < times["outside"]:
        raise AssertionError("order inside < middle < outside not reproduced")
    print("Same U-value for all three; insulation inside reaches 18 C %.0fx sooner than outside."
          % (times["outside"] / times["inside"]))
    print("EVERYTHING REPRODUCES")


if __name__ == "__main__":
    run()
