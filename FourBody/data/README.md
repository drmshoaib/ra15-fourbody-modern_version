# Initial condition files

All files use the format expected by `FourBody.exe`:

```
x0 y0 x1 y1 x2 y2 x3 y3        (positions)
vx0 vy0 vx1 vy1 vx2 vy2 vx3 vy3 (velocities)
m0 m1 m2 m3                      (masses)
T                                 (total integration time)
```

Units: G = 1, dimensionless.  No comment lines — values only.

---

## initialcond.txt — unit square (R = 1)

Four equal masses (m = 1) at the corners of a diamond of radius 1:

| Body | x | y |
|------|---|---|
| 0 | +1 | 0 |
| 1 |  0 | +1 |
| 2 | −1 | 0 |
| 3 |  0 | −1 |

Exact circular-orbit speed: v = sqrt(1/sqrt(2) + 1/4) ≈ **0.97832**

System maintains a rigid square for ~7 orbits (t ≈ 44) before chaotic breakdown.

---

## square_large.txt — larger square (R = 2, T = 200)

Same geometry scaled to radius 2.  Speed scales as v ∝ R^(−1/2):

v(R=2) = sqrt(0.9571 / 2) ≈ **0.69177**

Slower orbit; instability arrives later and trajectories spread over a wider region.

---

## square_tight.txt — tight square (R = 0.5, T = 30)

Compact version at radius 0.5:

v(R=0.5) = sqrt(0.9571 / 0.5) ≈ **1.38355**

Higher speeds and shorter inter-body distances accelerate instability onset.

---

## hierarchical.txt — hierarchical binary-binary (T = 300)

Two close binaries (pair separation = 1) whose centres of mass orbit each
other at a separation of 8 units.

```
Pair 1 (bodies 0, 1): centred at (+4, 0)
Pair 2 (bodies 2, 3): centred at (−4, 0)
```

Velocity derivation:

| Contribution | Formula | Value |
|---|---|---|
| Within-pair orbital speed | v_bin = sqrt(m / 2) = 1/sqrt(2) | 0.70711 |
| Pair-CoM orbital speed | v_out = sqrt(G * M_total / (4 * R)) = sqrt(4/32) | 0.35355 |

Each body's velocity = binary-pair velocity + pair-CoM drift velocity.

The integrator conserves energy to ~3 × 10⁻¹⁴ over T = 300 — a
near-integrable system well-suited for precision benchmarking.
