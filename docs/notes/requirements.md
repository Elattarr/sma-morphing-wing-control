# Requirements (draft v0.2) - standalone free wire
T = thesis, D = derived, A = my assumption (change freely, state it in docs/ASSUMPTIONS.md).

| ID | Requirement | Value | Source |
|----|-------------|-------|--------|
| R1 | Steady-state tracking error of contraction | <= 0.1 mm | A |
| R2 | Settling time (2 %) for a step to full contraction | <= 10 s | A (thesis heating cases 2.5-10 s) |
| R3 | Overshoot past the commanded contraction | 0 % (recovery is irreversible, one-way) | D |
| R4 | Maximum wire temperature | <= 80 C | A |
| R5 | Energy per full contraction (22 -> ~52 C) | report value from sim (thesis 22 -> 73 C: 173-199 J) | T |
| R6 | Peak current | <= 35 A until the real supply limit is known | T (placeholder) |
| R7 | Plant vs ANSYS recovery map | <= 0.25 K at equal recovery (now ~0.19 K) | D |
| R8 | Ambient range | 15-35 C | A |
| R9 | References | monotonic non-decreasing only (no return stroke) | D |
| R10 | Sensing | position noise 0.02 mm (1 sigma), sampling <= 100 Hz | A |
| R11 | Robustness | R1-R3 hold for h 5-82 W/m2K, c +-20 %, resistance +-20 %, ambient per R8 | A |

Full contraction: ANSYS 6.64 mm (3.32 % of 200 mm); frozen-strain estimate 6.54 mm (3.27 %). The plant uses 6.64.
Free-wire transitions (MEFF, T0 = 22 C): As = 50.79 C, Af = 51.82 C, window h*eps_L/beta = 1.04 K (see MEFF_window_findings.md).
Consequence: ~30 K of heating produce no output; ~98 J latent heat is absorbed inside a 1 K window.
Scope: heating only. No return stroke. No load on the wire during heating.
