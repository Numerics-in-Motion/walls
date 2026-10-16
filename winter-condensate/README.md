# Same brick, same insulation, same U-value: which wall collects the most condensate over the code's winter?

The three model walls of this repository's first study, unchanged (200 mm of solid brick, λ 0.77 W/mK, and 100 mm of
mineral wool, λ 0.035; insulation **inside**, **in the middle** or **outside** the brick; the same U-value), put
through the **DIN 4108-3:2001 Annex A** Glaser block calculation: steady vapour diffusion over a 1440 h winter
period, inside 20 °C / 50 %, outside −10 °C / 80 %, R_si 0.13, R_se 0.04 m²K/W. Endpoint: **water-equivalent
condensate** accumulated in the wall at the end of the winter period (kg/m²).

| build (outside → room) | condensate after 1440 h, brick μ 10 | over every brick μ 5–100 per leaf |
|---|---|---|
| insulation inside: brick 200 · wool 100 | **8.02 kg/m²** | 7.96 – 8.07 |
| in the middle: brick 100 · wool 100 · brick 100 | **0.67 kg/m²** | 0 – 1.39 |
| insulation outside: wool 100 · brick 200 | **0** | 0 (≥ 60.5 Pa below saturation everywhere) |

**Claimed:** the inside-insulated wall collects **at least 5.74×** as much as either other wall, at the joint
adverse corner (least inside value against the greatest middle value over every brick μ cell). **Not claimed:** an
order between the middle and the outside-insulated walls. With ordinary clay brick (μ 5, 10 or 16) the middle wall
always condenses (0.32–1.34 kg/m²) and the outside-insulated wall never does, but for some tested
combinations with clinker-type brick (μ 50/100) in the middle wall's inner leaf it also collects zero, so that pair is
left **unresolved**.

Why: with the insulation on the room side, vapour crosses the open mineral wool (μ 1) almost unhindered and meets
the brick's cold inner face, which lets it out only slowly (brick μ 10). With the insulation outside, the brick slows
the vapour while it is still warm and the open wool on the cold side lets it out — the line never reaches saturation.

```
python reproduce.py
```

`reproduce.py` re-runs everything below and compares every number with `reference/result.json`; it exits nonzero,
naming the quantity, if anything disagrees or a source changed (pure Python, under a second).

## Method

* Each layer is split into sub-layers of R ≤ 0.25 m²K/W (ISO 13788:2012 6.4.1, used here as an implementation
  choice), each refined ×16; the unrefined sub-layering is kept as a check (largest change 0.09 %).
* Temperatures follow the thermal resistances; saturation pressure from DIN 4108-3:2001 A.4.2
  (288.68 (1.098 + t/100)^8.02 at or above 0 °C, 4.689 (1.486 + t/100)^12.30 below).
* Vapour-pressure line = the lower convex hull of the saturation pressures in the (s_d, p) plane between the outdoor
  and indoor air values (the Glaser construction). Where it touches saturation, condensation occurs at a rate equal to
  the inflow minus the outflow, with δ₀ = 1/1.5·10⁶ kg/(m h Pa) (Z = 1.5·10⁶ s_d m² h Pa/kg).
* An independent envelope construction (minimum forward slopes) gives the same condensate in every cell (difference 0).
* Brick vapour-resistance factor μ is not part of the rig's thermal model, so it is swept: 5, 10, 16, 50, 100 for each
  brick leaf independently (DIN 4108-4 solid/perforated clay brick 5/10, ISO 10456 fired clay 16/10, DIN 4108-4 clinker
  50/100). Mineral wool μ 1.
* The pair rule: a pair is ordered only if its gap exceeds max(10 %, 3u), u = numerical uncertainty within the
  prescribed calculation (refinement 0.09 % + outside-check residual 0.006 % + saturation formula 0.005 %). The inside
  wall's gap is 474 %.

## Outside check

Xella/Hebel technical handbook, Example 3 (250 mm AAC μ 5 λ 0.14 + 9 mm bitumen μ 20000, R_si 0.25 / R_se 0.04,
−5/80 and 20/50 over 90 days): printed condensate **0.905** kg/m² and evaporation capacity **1.002** kg/m²; this solver
gives **0.9049** and **1.0022**.

## Controls

* Every layer μ 10 (mechanism diagnostic): inside 0.75, middle 0.33, outside 0 kg/m² — the order already follows from
  where the cold plane sits; the wool–brick μ contrast multiplies it.
* The outside-insulated wall with a vapour-tight outer skin (9 mm bitumen, μ 20000) condenses 0.41 kg/m², so its zero
  is not a blind spot of the method.
* A humid interior (20 °C / 80 %): inside 18.4, middle 1.29, outside 0 kg/m².
* A later block climate (the 2014 edition's: 90 days, −5 °C / 80 % outside, R_si 0.25, δ₀ = 2·10⁻¹⁰ kg/(m s Pa)) gives
  the same order in every cell (inside 10.4–10.6 kg/m²). This is a **boundary-condition sensitivity check only**: the
  2014 and 2018 editions restrict the simplified method for interior insulation with R > 1 m²K/W, so no statement is
  made here about which edition applies to which wall.

## Scope

* A code's simplified block calculation (steady diffusion, no capillary transport, no sorption, no rain, no air
  movement) applied to this rig's idealised walls; the condensate is a water-equivalent amount at planes below 0 °C.
  It is not a prediction of moisture in a real wall, not a verdict on any construction, and drying in the following
  season is not compared.
* The scene is the 2001 edition, which has since been superseded.

## Sources

* DIN 4108-3:2001 provisions and Annex A (block climate, R_si/R_se, Z, saturation formula) — umwelt-online,
  a third-party HTML transcription:
  <https://www.umwelt-online.de/regelwerk/bau/din/4108_31.htm>,
  <https://www.umwelt-online.de/regelwerk/bau/din/4108_32.htm>,
  <https://www.umwelt-online.de/regelwerk/bau/din/4108_33.htm>
* Outside check and δ₀ for the later block: Xella/Hebel technical handbook, Example 3 —
  <https://technik.xella.de/media/books/handbuch/151/>, <https://technik.xella.de/media/books/handbuch/152/>
* Vapour-resistance factors: DIN 4108-4 table (secondary transcription)
  <https://www.scribd.com/document/767327809/DIN-4108-V-4>; ISO 10456:2007 table (secondary transcription)
  <https://www.scribd.com/document/736512244/10456-2007-ISO>
* Later editions' restriction for interior insulation: DIN Media, DIN 4108-3:2014 —
  <https://www.dinmedia.de/de/norm/din-4108-3/222022641>; later block values (secondary):
  <https://de.wikipedia.org/wiki/Glaser-Verfahren>
* The walls: this repository's first study (`../reference/result.json`).
