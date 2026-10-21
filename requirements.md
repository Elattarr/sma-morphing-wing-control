# Requirements (draft v0.1) - standalone free wire
T = taken from thesis, A = my assumption (change freely, state it in the README).

| ID | Requirement | Value | Source |
|----|-------------|-------|--------|
| R1 | Steady-state tracking error of wire contraction | <= 0.1 mm | A |
| R2 | Settling time (2 %) for a step to full contraction | <= 10 s | A (thesis heating cases 2.5-10 s) |
| R3 | Overshoot | <= 5 % of step | A |
| R4 | Maximum wire temperature | <= 80 C | A |
| R5 | Energy per full contraction (22 -> ~52 C) | report value | T (thesis energies for 22 -> 73 C: 173-199 J) |
| R6 | Peak current | <= 35 A | T (2.5 s case) |
| R7 | Output model error vs ANSYS | <= 5 % of full contraction | A |

Full contraction target: 3.27 % x 200 mm = 6.54 mm (frozen strain from the thesis).
Free-wire transitions (derived from MEFF constants, T0 = 22 C): As = 50.8 C, Af = 51.8 C (window h*eps_L/beta = 1.04 K, see MEFF_window_findings.md).
Scope: heating only. No return stroke. No load on the wire during heating.
