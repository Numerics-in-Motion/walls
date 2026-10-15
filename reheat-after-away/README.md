# Back from a night away, or a month away: which room is warm first?

The three model rooms of this repository's first study (same brick, same insulation, same U-value; insulation
**inside**, **in the middle** or **outside** the brick) are held at 20 °C with 0 °C outdoors. The heating goes off
while you are away; on return the same **1.5 kW plug-in heater** goes on. Endpoint: minutes until the room **air**
first returns to 20 °C.

| away | insulation inside | in the middle | insulation outside |
|---|---|---|---|
| 12 hours | 15.6 min | 2.0 min | 1.3 min |
| 28 days | 16.5 min | 22.7 h | 38.2 h |

**Claimed, pair by pair** (a pair's order is claimed only if its gap exceeds max(10 %, 3u) at its own worst-case
input-sensitivity index u, which includes furniture heat capacity 0–180 kJ/K as a structural assumption):

* **outside vs inside swap places**: after 12 h the outside-insulated room's air returns to 20 °C before the
  inside-insulated room's (the inside room takes ×12 as long, empty room); after 28 days the inside-insulated room
  returns before the outside-insulated one (the outside room takes ×139 as long). Gap vs need: 1141 % vs 776 % (12 h), 13817 % vs 457 % (28 d).
  The order flips in none of the 1,024 shared lo/hi input corners, with or without furniture.
* after 12 h the middle room is **within uncertainty** (not ranked); after 28 days all three rank inside < middle < outside.
* they swap after **about 1–4 days away** — between 22 and 93 h over the tested shared input corners with the
  1.5 kW heater fixed (nominal 42 h). This is a tested envelope, not a confidence interval or a proven bound.

**Depends on the heater** (swept 0.4–3 kW): after 12 h the inside-insulated room leads below about 496 W (no
reversal); between about 496 W and 903 W the pair is within uncertainty (a heater sized to the house, 657 W =
2× the steady loss, gives 57.0 vs 59.3 min, a tie); from about 903 W to 3 kW the outside-insulated room leads
and clears. 1.5 kW is one of the common ratings of portable electric heaters (sources in `reference/result.json`,
`heater_basis`).

The headline ratios are for the empty room; with 180 kJ/K of furniture on the air node they are ×7.6 (12 h) and
×35 (28 days), same order.

```
python reproduce.py
```

## Method
The parent study's finite-volume walls (`../reproduce.py`, unchanged; its frozen reference is reproduced first).
Initial state: the steady profile at 20 °C held by the steady loss; heater 0 W during the absence; then 1.5 kW.
Recovery times by an exact modal solution of the same finite-volume network, checked against the rig's
Crank–Nicolson stepping (≤ 1e-4) and a Laplace-domain superposition of the rig's own transfer-matrix solution
(≤ 1e-3); with massless walls all three rooms recover in the same time (control).

## Scope
An ideal model: six identical one-dimensional faces, no windows, leakage, ground or sun; one well-mixed air node;
air temperature only — not comfort, energy or cost. **Not retrofit advice**: interior insulation carries moisture
and freeze–thaw risks this dry model does not assess.
