# Bench protocol and hardware specification: current-driven one-way SMA wire

Status: everything in `docs/eval_*.txt` is **verification** (controller vs. the MEFF model it was built from, plus structural-mismatch plants).
Nothing has been **validated** on a real wire. This document is the path from the first to the second.

Reference configuration (what the evaluation assumed): 0.5 mm wire, 200 mm free length, bias spring 1.25 N/mm, preload at the extended end
290 MPa (57 N), current limit 3 A, 20 Hz control, ambient 22 °C.

## 1. Hardware specification

| Item | Requirement | Why (number from the evaluation) |
|---|---|---|
| Wire | 0.5 mm, MEFF constants of specimen C; 2 mm is impractical (needs 660–790 N preload to return at ambient) | `spring.py`: return at ambient needs σ_ext ≥ 211 MPa |
| Preload / spring | σ_ext ≥ 290 MPa at the extended end; spring rate ≈ 1.25 N/mm (measure it) | at 250 MPa a spread ≥ 3 K wire never returns fully; homing refuses (`preload_too_low`); at 290 MPa spreads up to 3 K (and 25 % tail to 20 K) pass, a 5 K spread is still refused |
| Current source | 0–3 A, ≥ 12 bit, hardware limit at 3 A independent of software | 3 A is [ASSUMED] for 0.5 mm; contraction in ≈ 6 s needs ≈ 2.7 A |
| Current / voltage measurement | I ±0.5 %, V at the wire (4-wire), both at ≥ 20 Hz | monitor checks resistance and current mismatch; identification needs R to 1 % |
| Control thermocouple T1 | type K, ≤ 0.1 mm, bonded at mid-length, lag ≤ 1 s, noise ≤ 0.3 K | evaluated up to 0.5 K noise and 3 s lag (4/4 targets) |
| Independent channel T2 | second thermocouple at a different position, own amplifier, read by the **safety** processor | single-sensor monitoring does NOT catch a frozen thermocouple (see `eval_thin_safety.txt`) |
| Position sensor | noise ≤ 0.02 mm rms (0.1 mm tolerated), ≥ 20 Hz, range ≥ 10 mm (stroke = 3.3 % of 200 mm ≈ 6.6 mm) | evaluated up to 0.1 mm noise and 0.05 mm quantisation |
| Safety relay | normally-open, in series with the wire, driven by the safety processor only | trip = relay open = power off; wire then extends under the spring (safe end) |
| Safety processor | separate MCU or comparator logic, own supply, receives I_cmd, I_meas, T1, T2, d, V and a controller heartbeat | `safety.py` runs on it; controller failure must not disable it |
| Thermal fuse (recommended) | cut-out near the wire clamp, ≈ 130 °C | last line of defence if both channels and the processor fail |

## 2. Commissioning gates (do not skip, do not reorder)

**Gate 0, electronics, no wire.** 1 Ω power resistor in place of the wire. Verify current loop, I and V readings, the 3 A hardware limit, relay opening.
Fault injection on the dummy load: command 0 A while the stage is forced on (expect `current_mismatch` within 1 s and relay open); freeze the controller
heartbeat (expect `watchdog` within 0.5 s); disconnect each thermocouple (expect `sensor_range` or `tc_disagree`).

**Gate 1, sensors.** Thermocouples against a reference thermometer at ambient and ≈ 60 °C (±1 K). Position sensor against gauge blocks. Record 60 s at rest:
the noise values set the homing motion threshold (done automatically) and tell you which row of `eval_*_sensors.txt` applies.

**Gate 2, thermal identification (Test A).** Wire clamped with the spring, extended, cold. 10 s at rest, constant current of about 1 K/s initial rate,
cut when T1 rises 12 K (stays far below A_s), cool to within 1.5 K of ambient. Run `fit_thermal` (`identify.py`).
Accept: residual rms ≤ 2× thermocouple noise; m·c within ±30 % of 0.21 J/K, h·S within a factor of 2 of 9.4 mW/K (otherwise find out why before going on).
The fitted numbers go into `SMAParams` (c, h) and `SafetyCfg.h_min/h_max` is checked against the fitted h·S.

**Gate 3, homing at reduced power.** First run with I_max = 1.5 A and T_abs_max = 110 °C. Then 3 A and 125 °C.
Accept: no fault; A_s, A_f, T_fs, T_fe plausible for the stress (MEFF: A_s 87 °C, A_f 94 °C, T_fs 35 °C, T_fe 28 °C at 250 MPa); measured stroke ≥ 30 % of 6.6 mm;
preload margin T_fe − T_amb ≥ 3 K (the code refuses below 2 K; take 3 K as the bench rule because the ambient drifts).
Then run Test B/C (`fit_latent`, `fit_resistance_contracted`) on the same record: m·Λ within ±40 % of 6.1 J, R_A/R_M ≈ 1.25.

**Gate 4, monitor on the real wire.** With the wire connected and power limited to 1.5 A: repeat the Gate 0 injections (T1 unplugged, T2 unplugged,
position sensor unplugged during a heating step, controller stopped). Every injection must open the relay and leave the wire below 130 °C.

**Gate 5, closed loop.** Small steps first (r = 0 → 0.25 → 0), then the evaluation sequence r = 1.0 → 0.5 → 0.8 → 0.0 with a 30 s hold at each.
Compare with the evaluation: hold error ≤ 2 %, no overshoot beyond 5 %, no monitor trip, peak temperature ≤ 125 °C.
Anything outside that is a model error: log it, do not retune around it (see §4).

**Gate 6, repetition.** 100 cycles. Log A_s, A_f, stroke at every re-homing. MEFF has no functional fatigue; a real wire does, so a drift here sets the re-homing interval.

## 3. Logging for the thesis (validation data)

20 Hz, one CSV per run: t, I_cmd, I_meas, V, T1, T2, d, controller mode, monitor state. From the homing runs extract A_s, A_f, T_fs, T_fe and the strains and compare them with the MEFF
prediction at the measured stress; this comparison (model vs. real wire) is the validation the evaluation could not provide.

## 4. Assumptions the bench must settle (all [ASSUMED], none measured)

- Convection coefficient h = 30 W/m²K, extrapolated; identification (Gate 2) replaces it. Evaluated range 5–82.
- Current limit 3 A and T_abs_max = 125 °C for the 0.5 mm wire.
- Thermocouple lag ≤ 4 s (monitor `tc_lag_max`), h bounds 5–100 W/m²K, ambient ± 10 K (monitor).
- Uniform temperature along the wire (lumped model). A real wire has end effects at the clamps.
- Ambient 15–26 °C and air-flow factors 0.5–3 were tested; outside that nothing is known.
- Window shape: the plants tested were single-element and parallel-element (spread up to 5 K, tail to 20 K). A wire that does not behave like a sum of such elements is untested.
