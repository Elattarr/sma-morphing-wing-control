# Status: standalone one-way SMA wire controller, simulation stage

Date of this status: 2026-10-30. Scope: free one-way wire pulled back by a bias spring, driven by current. Python reference, portable C99 port, simulation only.

## Done (verified against the MEFF model and structural-mismatch plants, not against a real wire)

| Item | Where | Evidence |
|---|---|---|
| Plant: MEFF ratchet + energy balance + bias spring (coupled stress), parallel-element plant for wide windows / tails / sub-loops | `src/sma_control/{cycle_plant,spring,elements}.py` | `tests/test_spring.py`, `tests/test_elements.py` |
| Bidirectional position controller (hold / up / down), overshoot exit, adaptive limit and hold temperature | `cycle_controller.py` | `docs/eval_*` |
| Homing / calibration (power-up cycle): measures stroke end points and A_s, A_f, T_fs, T_fe; corrects wrong stress, wide windows; refuses when the preload is too low | `HomingRoutine` | `tests/test_homing.py`, robust and struct tables |
| Independent safety monitor: hard over-temperature (two channels), channel disagreement, energy-balance envelope, stall, current mismatch, wire resistance, watchdog, latching relay | `safety.py` | `tests/test_safety.py`, `docs/eval_*_safety.txt` (fault injection, with/without monitor), zero false trips in every other run |
| Parameter identification from bench data (m c, h S, thermocouple lag, R_M, R_A, latent heat) | `identify.py`, `scripts/identify_bench.py` | `docs/identify_bench.txt`, `tests/test_identify.py` |
| Bench protocol and hardware specification with go/no-go gates | `docs/BENCH_PROTOCOL.md` | - |
| C99 port of controller, homing and monitor; replay of logged Python runs sample by sample (double: identical to 4e-14 A; single: same trips and calibration, Tfe within 0.15 K / 1.5 K noisy) | `c/` | `tests/test_c_port.py` |
| Other alloys with unchanged gains (12 MEFF variants, nominal and own constants) | `scripts/evaluate_materials.py` | `docs/eval_thin_materials.txt`, `tests/test_materials.py` |

## Not done (needs the real wire or is outside this stage)

1. Validation on hardware: every number above compares the controller with the model it was built from (or with a perturbed copy of it).
2. The wing: coupling to the airfoil stiffness, aerodynamic load, the embedded wire map.
3. The two-way effect (not in MEFF here).
4. MATLAB/Simulink port.
5. Real-time timing, ADC/DAC, EMC and power-stage behaviour of the C code on a target MCU (the code is portable C99, compiled and replayed on a PC only).

## Known limits found in the evaluation

- A single thermocouple that freezes is **not** caught by the envelope check alone on the 0.5 mm build wire (2 of 2 runs not stopped, Tmax 116.7 C, `eval_thin_safety.txt`); it is caught on the 2 mm wires. With the independent second channel all 14 two-channel fault runs stop. T2 is mandatory.
- Thermocouple lag is identified poorly when it is small (0.1 s true -> 0.047 s fitted, 0.5 s -> 0.45 s); m c, h S, R, m Lambda are within 1 %. The lag only enters the monitor through a fixed bound (4 s), so this does not matter for safety.
- The preload must cover the transformation-temperature spread of the wire: at 250 MPa a spread of 3 K or more (and a 12 K window) never returns fully and homing refuses; 290 MPa covers spreads up to 3 K and the tail cases; a 5 K spread (window about 16 K wider) is still refused (`eval_thin_struct_preload290.txt`). A wire with that spread needs a higher preload, which also raises the peak temperature (122 C already reached of 125 C).
- Ambient 26 C at 250 MPa leaves 1.7 K return margin and is refused; at 290 MPa it works.
- Homing costs 200-900 s per power-up (0.5 mm wire; up to 3400 s for a 2 mm wire at h = 5) and heats the wire to A_f + ~17 K; peak temperature reaches 122 C of the 125 C limit for some alloys and wide windows.
- Single-precision C: threshold decisions can shift by one sample (accumulated `still += Ts`); tolerable, but use integer sample counters if bit-level agreement is wanted.
- The 2 mm wire cases are kept for the thesis wire only; the build recommendation is 0.5 mm.

## Reproduce

```
python -m pytest                                    # all tests
SMA_D_MM=0.5 SMA_IMAX=3 SMA_SPRING_K=1.25 SMA_SIG=290 python scripts/evaluate_controller.py robust|struct|sensors|steps|ramps|disturb|safety
python scripts/evaluate_materials.py ; python scripts/identify_bench.py
python scripts/export_replay.py out/ && make -C c && c/replay out/nominal.txt
```

## Embedded-wire plant (first run, 2026-10-30)

Plant: 2 mm wire, 200 mm, bias spring 5 N/mm + wing structure about 6.4 N/mm at the wire (derived from the thesis ANSYS Up 40 mm run, force 7 N -> 48 N for 3.27 % recovery; assumes the spring acts in parallel with the structure at the wire end), preload 250 MPa (785 N), current limit 35 A, ambient convection h = 30 W/m2K (the real heat loss into the CFRP is unknown). Results: `docs/eval_wing_*.txt`, figures `docs/figs_wing/`.

- Windows with the spring: A_s 87.3, A_f 91.8, T_fs 32.2, T_fe 27.7 C (thesis ANSYS without spring: about 51-55 C).
- Contraction 7 s to full stroke at 35 A, peak 92 C, 350 J; extension 90-380 s (cooling limited).
- Sequence 1.0 -> 0.5 -> 0.8 -> 0.0: 4 of 4, hold error 0.8 %, about 1000 s, homing 690 s, peak 110 C.
- Fails or refuses: ambient 26 C (and +8 K during a move), 12 K window (`preload_too_low`); h = 5 reaches 3 of 4 in 5000 s.
- Not validated: the structure stiffness comes from one ANSYS run; the thermal model is air convection, not CFRP.

## Heat-only (one-way) controller on the thesis wire (2026-10-30)

Plant: 2 mm × 200 mm NiTi, structure only (6.4 N/mm at the wire, 7 N = 2.2 MPa preload), no bias spring, 35 A limit, no homing
(homing would use the one-way stroke). TE map: from the ANSYS Up 40 run (see 'ANSYS comparison' below), not linear. Code: `cfg_heat_only` in
`cycle_controller.py` (no 'down' mode, hold = power off, pre-heat to A_s − 3 K, 1 K/s temperature creep until motion is seen, then
position loop Kr 8/s, reference rate 0.1/s); harness `scripts/evaluate_heatonly.py`; results `docs/eval_heatonly_*.txt`, `docs/figs_wing/heatonly_*.png`.
All of it is model vs model (verification), nothing validated.

What the runs show (default preset, regenerated 2026-10-30 after the last fixes; older numbers such as "14 of 17" are superseded):
- Steps from r = 0: full 42.2 mm settles in 14.5 s, overshoot ≤ 0.7 %, steady-state error ≤ 0.3 %, peak wire temperature 54.2 °C, about 220 J per full stroke.
- Staircase 0.25 → 1.0: 4/4, hold error 0.2 % nominal. Disturbances (air ×3, ×0.5, ambient ±10 K) change the final position by ≤ 0.4 %.
- Ramps 0 → 1: RMS error 6.7 % (20 s), 2.4 % (60 s), 1.8 % (180 s and 600 s); max error 26 % on the 20 s ramp (the loop is deliberately slow: reference rate 0.1/s).
- Parameter errors: 14 of 16 cases reach 4/4 (worst hold error 0.83 %), including the plant fitted to the ANSYS Up 40 data (hold error 0.70 %, Tmax 56.5 °C). Two reach 3/4: window spread 5 K, and spread 2 K + 25 % tail to 20 K. Both are far wider than the ANSYS data (spread 0.4 K, 5 % tail over 5 K), so they are stress tests, not the expected case.
- Sensors: 6 of 7 pass (noise 0.5 K, lag 3 s, 0.05 mm quantisation, "all" combined). Displacement noise 0.1 mm fails (0 of 4). Measured on a cold wire:
  false motion detection at 0.02 mm noise 0/10 seeds, 0.03 mm 1/20, 0.05 mm 4/20. Requirement for the bench: displacement noise ≲ 0.02–0.03 mm on the wire stroke
  (about 0.4 % of 6.5 mm); LVDT/laser on the wire end, not a TE camera. Motion detection defaults (creep_confirm 3, fast_detect 3.0) were chosen by a sweep over the
  robust + sensor sets; they are the least invasive variant, not an optimum.
- Safety fault injection (30 runs): all end below trip + 10 K; trips: watchdog, current_mismatch, resistance, tc_disagree, tc_reads_low, stall (frozen position sensor at 28 s).
  Frozen thermocouple (bit-identical T1 for 10 s while heated or hot) is now caught as `tc_frozen` in all 3 injection times (30 of 30 runs end safe, every injected fault trips); zero false trips in the other
  heat-only sections and in the 112-test suite. A frozen reading that still carries noise-like variation, or a drifting sensor, is not covered.
- Found and fixed while porting: the stall check counted a hot wire cooling at rest as "no motion" and tripped after the learner stop (wire overshoots to ~65 °C, then cools ≥ 8 s above
  A_f,max + 5 K). The stall timer now only counts while I_meas > 0.3 A (Python and C).

Controller features added after the first version (all model-only): stop at the end of travel (no creeping to the temperature cap; sticky accepted target), one-point zero at
power-up (assumes the wire is extended at start), first-stroke learner (two-point calibration from the measured end positions, A_s from the first motion, A_f from the
temperature at the last motion), faults `position_dropped`, `stalled_short`, `reset_timeout`, `begin_reset()` re-arm, `rate_dn = 0` (a negative rate_dn had kept the
reference above the target, so "done" never fired; this was the biggest single improvement).

Support code: `te_map.py` (r ↔ TE mm, piecewise linear, default 0 → −42.2 mm, ASSUMED linear), `hardware_io.py` (`RigIO` interface + `run_step`: sensors → controller →
safety monitor → permitted current and relay; untested on hardware), C port (`c/sma.c`, `c/replay.c`): bit-identical to Python in double precision on six heat-only replays
(nominal, 4 K spread, noisy, 10 %-short wire with learner, two fault injections); float build agrees on decisions, current differs by ≤ 3e-3 A.
Tests: `tests/test_heat_only_unit.py` (12 behaviour tests), `tests/test_heat_only.py`, heat-only cases in `tests/test_c_port.py`; full suite 112 passed, 1 xfailed.

ANSYS comparison (2026-10-30, run heating 22 → 80 °C, `docs/ansys_up40_recovery.csv`; recovery r from the wire-end displacement, stroke 6.48 mm; model stroke 6.48 mm):

| r | model T [°C] | ANSYS TE [°C] | ANSYS wire end [°C] |
|---|---|---|---|
| 0.1 | 51.5 | 51.5 | 51.3 |
| 0.5 | 52.7 | 52.7 | 52.4 |
| 0.9 | 53.9 | 54.3 | 54.1 |
| 0.99 | 54.2 | 55.9 | 55.9 |
| 0.999 | 54.2 | 59.6 | 60.8 |

- The model matches the 10–50 % points within 0.3 K and the stroke within 0.1 %; its window ends sharply at 54.2 °C, ANSYS has a tail (99 % at 55.9 °C, 99.9 % at 60–61 °C). The plant used in the evaluation (parallel elements, spread 0.4 K, 5 % tail over 5 K, RMS 1.8 % of the stroke vs ANSYS) reproduces that tail; the controller reaches 4/4 on it.
- The earlier "recovery to 73–77 °C" is the local equivalent-plastic-strain MAXIMUM (hot spot): 0.033 at 50 °C, 0.003 at 59 °C, 0.0003 at 71 °C, 0 at 76 °C. The whole-wire recovery (wire end and TE deflection) is 99.9 % at 60 °C. The local max total strain stays at 0.35 % at 80 °C (residual).
- TE deflection −6.63 mm (preloaded) → −42.24 mm; the TE curve lags the wire end by up to 13 % of the stroke (concave): the linear map was wrong in the middle. `te_map.py` now uses the ANSYS points.
- The ANSYS run has no average-strain output for the wire in this file (only min/max); the end stress is the local max (178 MPa at the end), not the 15 MPa mean I assumed earlier. Not reconciled.

Still open / not validated:
- Everything above is model vs model or model vs ANSYS (a simulation); no bench data. Sensor noise requirement ≲ 0.03 mm, gains and safety thresholds are unvalidated.
- Supervisor A_s/A_f 51/69 °C is not consistent with the ANSYS whole-wire window (51–56 °C, tail to 60 °C); thermal model is air convection, not CFRP; resin T_g vs 56–60 °C (up to 76 °C locally) is unchecked.
- Return mechanism is not in the controller (power off, wait for the displacement sensor, timeout, `begin_reset`); the wing/aero load that would extend the wire is not modelled.
- The TE map is from one load case; the one-point zero (wire extended at start) is an assumption.
