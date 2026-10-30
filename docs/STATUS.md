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
