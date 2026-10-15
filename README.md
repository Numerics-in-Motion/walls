# Same wall, same U-value: which room warms first?

Three model rooms are built from the same 200 mm of solid brick and 100 mm of insulation per
square metre of envelope. Only the position of the insulation changes, so all three have the
**same U-value (0.3042 W/m²K)** and the same steady heat loss. A 1.5 kW heater
switches on in a room that starts at 0 °C throughout. How long until the room **air** reaches 18 °C?

| build (outside → room) | air reaches 18 °C | with 180 kJ/K added to the air node |
|---|---|---|
| inside: brick 200 · insulation 100 | **14.1 min** | 58.5 min |
| middle: brick 100 · insulation 100 · brick 100 | **19.52 h** | 20.31 h |
| outside: insulation 100 · brick 200 | **34.52 h** | 35.34 h |

Equal steady heat loss does not mean equal switch-on time: with the insulation on the room side,
the brick stays out of the warm-up; with it outside, the room-side brick absorbs much of the
heater's energy as the air warms.

```sh
python -m pip install -r requirements.txt
python reproduce.py
```

`reproduce.py` solves each wall with cell-centred finite volumes (Crank–Nicolson, 1.25 mm cells,
1 s steps), closes the energy books, checks every finish time against an independent
Laplace-domain solution (layer transfer matrices, Talbot inversion), and compares each number with
`reference/result.json`. It exits nonzero, naming the quantity, if anything disagrees.

## The model (an ideal comparison, not a representative cabin)

* A 3 m interior cube whose six faces all carry the same one-dimensional wall; no windows, air
  leakage, ground contact, corners, sun or moisture. Outdoors fixed at 0 °C.
* Brick λ 0.77 W/mK, ρ 1700 kg/m³, c 800 J/kgK; insulation λ 0.035 W/mK, ρ 25 kg/m³,
  c 1400 J/kgK — representative values.
* Surface resistances R_si = 0.13 and R_se = 0.04 m²K/W on every face, used as fixed model inputs
  (the wall convention; roof and floor values differ).
* Room air is one well-mixed node; the result is an **air temperature**, not comfort.
* The heater is constant at 1.5 kW. The steady loss at 20 °C in this model is only about
  329 W, so it is oversized for holding a temperature, as plug-in heaters are.

## What this is not

* **Not retrofit or building advice.** Interior insulation carries moisture and freeze–thaw risks
  that this dry-wall model does not assess. No comfort or energy-saving claim is made.
* A split build (50 mm of insulation on each side of the brick) was also run; the tested 50 mm
  room-side layer made it effectively tie the inside build, so it is not one of the three.

## More studies on this rig

* [`reheat-after-away/`](reheat-after-away/) — the same three rooms after a night away and after a month away: which is warm first?
