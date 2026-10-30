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
(homing would use the one-way stroke). TE = r × 42.2 mm is an ASSUMED linear map (thesis Up 40 mm). Code: `cfg_heat_only` in
`cycle_controller.py` (no 'down' mode, hold = power off, pre-heat to A_s − 3 K, 1 K/s temperature creep until motion is seen, then
position loop Kr 8/s, reference rate 0.1/s); harness `scripts/evaluate_heatonly.py`; results `docs/eval_heatonly_*.txt`, `docs/figs_wing/heatonly_*.png`.
All of it is model vs model (verification), nothing validated.

What the runs show (nominal numbers, clean sensors): full 42.2 mm in about 14 s settle, overshoot ≤ 2 % of stroke on steps,
staircase 0.25 → 1.0 reached 4/4 with hold error ≤ 0.6 % (0.25 mm TE), peak wire temperature 54.3 °C, about 220 J per full stroke,
disturbances (air flow ×3, ×0.5, ambient ±10 K) do not change the final position (≤ 0.8 %). Slow ramps track with RMS 2.1–2.5 % (about 1 mm TE).

Limits found (do not hide them):
- The first version stalled: after a small overshoot the old logic entered 'down', which does nothing in a one-way plant (deadlock). Fixed by `heat_only`.
- An overshoot cannot be undone. The result is very sensitive to how the wire enters the window: with the old pre-heat (integrator still heating at ~20 K/s)
  most parameter-error cases overshot by about 6 % (2.6 mm TE) at the first target. The creep phase removed that.
- Parameter errors (Λ, c, h, ρ, ambient, controller stiffness ×0.5/×2, preload error 2 → 30 MPa in the controller): 14 of 17 cases reach all 4 targets,
  worst hold error 1.6 %. Still failing: true preload 30 MPa vs 2.2 MPa assumed (first overshoot 2.4 %, just outside the 2 % band), and spreads of 5 K or a 25 % tail
  to 20 K (the last target needs more than the temperature cap allows; 70–75 °C, ~3 kJ).
- A tuning change that looked better (creep_detect 0.03) made 6 of the parameter cases worse; the controller is tuned at the edge. Treat the gains as a starting point for the bench.
- Displacement sensor: noise of 0.1 mm (1.5 % of the wire's 6.5 mm stroke) or 0.05 mm quantisation makes it fail (0 of 4). Needs noise and resolution
  below about 0.02 mm on the wire stroke (an LVDT or laser on the wire/end block, not a TE camera). Thermocouple lag 3 s and noise 0.5 K are tolerated.
- Return/reset is not in the controller (power off, wait for the displacement sensor, timeout, re-home) and the CFRP surroundings are not modelled.
