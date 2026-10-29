# Plant assumptions and the measurement that replaces each

| Parameter | Value | Source / tag | Verified? | Replace with |
|---|---|---|---|---|
| MEFF h, R, beta, eps_L, T0 | 145, 167, 5.6, 0.04, 22 | THESIS (Alazzawi specimen C) | yes | - |
| Lode m | 0.01 | user-supplied ANSYS input | yes | - |
| A_s, A_f, window | 50.867, 51.903 C, 1.036 K | MEFF with m = 0.01 (m = 0: 50.79 / 51.82) | yes (ANSYS recovery agrees within data resolution) | real-wire DSC (see window row below) |
| xi0, full stroke | 0.954, 6.64 mm | user / ANSYS max node | yes | - |
| d, L | 2 mm, 200 mm | user | yes | - |
| rho (density) | 6450 kg/m3 | Dynalloy datasheet, Alazzawi | yes | - |
| c | 837 J/kg/K | Dynalloy datasheet (0.2 cal/g/C); old value 320 came from an ANSYS example, unsupported | datasheet only | DSC / bench heating curve. Sweep 320-837 |
| latent heat | 0.156 J/mm3 = 98 J | Dynalloy datasheet 24.2 J/g, Shayanfard | datasheet only | DSC. Alt 46 J (Clausius-Clapeyron, T*beta*eps_L): swept |
| h_c | 30 W/m2K | ASSUMED, extrapolated to d = 2 mm from Pathak 2010 (68-153 for thin wires) and Battaglia 2025 (81.9 at 0.5 mm) | NO | bench current step, cooling tail. Sweep 5, 20, 82 |
| rho_e (M / A) | 0.8 / 1.0 uOhm m | Dynalloy datasheet | direction disputed (Pathak, Shayanfard: A < M) | 4-wire R cold and hot. Reversed direction swept |
| T_inf | 22 C | design | - | sweep 15-35 |
| noise_mm, quant_mm | 0.02, 0 mm | ASSUMED | NO | sensor datasheet / measured noise at rest |
| tail_fraction, tail_width | 0, 40 K | ASSUMED | NO | strain-vs-temperature curve of the real wire at ~0 stress |
| T_shift | 0 K | ASSUMED | NO | same curve |
| Imax | 35 A | THESIS R6 | - | actual supply/driver limit |

Unresolved: ANSYS forward-transformation stress rise (2.6 vs 7.1 MPa predicted) - not used by this plant.

## Robustness sweep (scripts/param_sweep.py; controller = nominal model, plant = perturbed)

Result file: the table printed by the script. Summary of what it showed:

* c 320-837, h_c 20-82, rho_e +-20 %, reversed rho_e, ambient 15-35 C: final error <= 0.015 (1.5 % of stroke), no overshoot > 4.4 %, Tmax <= 55 C.
* h_c = 5: overshoot 8.2 % at r_ref 0.5, Tmax 55.0 C at r_ref 1.0 (the R4 limit of 56 C is not violated).
* Lambda = 46 J instead of 98 J: overshoot 13.2 % at r_ref 0.5 (the controller assumes 98 J). Needs lam_scale tuning or latent-heat estimate.
* Window 12 K (same A_f): r_ref 0.5 overshoots to 0.97. The controller parks at A_s - 0.5 K of the nominal map, but this wire has already transformed there. One-way recovery cannot be reduced, so the overshoot is permanent. Needs the real map (bench) or a slow position-feedback approach.
* Known open item: tests/test_controller.py h_scale=0.2 case is xfail (3.2 % overshoot vs 3 % limit) since h_c default dropped from 81.9 to 30. Controller tuning is deferred by the user.

## Loading model (src/sma_control/loading.py, scripts/loading_study.py)

1D uniaxial MEFF with stress input. Compared with the ANSYS standalone run (T = 22 C, force ramp to 857 N, unload).

* Bulk transformation matches: 50 % of the wire transformed at 207.5 MPa average stress in both theory and ANSYS; 80 % at 209.7 (theory) vs 210.4 (ANSYS).
* Local onset at the clamp: s_max = 204.1 MPa (theory 203.97), while the average stress is only 164 MPa there: stress concentration 1.24 at the clamped end. A 1D model cannot reproduce this early local transformation (xi_avg = 0.05 at 189.5 MPa in ANSYS vs 204.3 in theory).
* xi_avg does not reach 1 in ANSYS (0.94 at 275 MPa): clamp regions with stress below the onset. The bulk is fully transformed.
* Free-end displacement: 1D is 2-3 % below ANSYS during loading and also 3 % below in the purely elastic range, so this is a constant node/geometry offset, not a MEFF effect. After unloading: ANSYS 6.629 mm, 1D (xi = 1) 6.532 mm.
* y_max = 6.64 mm therefore corresponds to a fully transformed bulk (xi about 1), not to the volume-average 0.954 / 0.940. Plant value unchanged.
* ANSYS has no thermal coupling: self-heating during fast loading is not in the data and is not validated here.

## Full-cycle controller, safety layer, identification, C port (added with the bias-spring work)

All numbers below are design choices or assumptions; none was measured on a real wire. Evaluation files: `docs/eval_*.txt`.

| Item | Value | Tag | Replace with |
|---|---|---|---|
| Wire for the real build | 0.5 mm, 200 mm; 2 mm needs 660-790 N preload to return at ambient | DESIGN | - |
| Bias spring | 1.25 N/mm (0.5 mm wire), 20 N/mm (2 mm wire) | DESIGN | measured spring rate |
| Preload (stress at the extended end) | 290 MPa = 57 N (0.5 mm). Rule: T_fwd(1) = T_ambient + 12 K, i.e. sigma = ((T_amb + 12 - T0) beta + R/(1+kappa m) + h eps_L)/sqrt(2/3) | DESIGN | load cell |
| Homing refusal | `preload_too_low` if the measured end of extension is < 2 K above ambient (the wire cannot return fully; the position scale would be wrong) | DESIGN | - |
| Current limit | 3 A (0.5 mm), 35 A (2 mm) | ASSUMED | actual driver / wire datasheet |
| T_abs_max | 125 C; homing peak is A_f + ~17 K (111 C nominal, up to 122 C for wide windows / other alloys) | ASSUMED | wire and insulation datasheet |
| Convection h | 30 W/m2K, evaluated 5-82 | ASSUMED | `identify.fit_thermal` (bench Test A) |
| Sensors (evaluation) | T1 noise 0.1 K (to 0.5 K), lag 0-3 s; displacement 0.02 mm noise (to 0.1 mm), quantisation to 0.05 mm; 20 Hz | ASSUMED | bench noise at rest |
| Safety monitor bounds | h in [5, 100] W/m2K, m c in [0.3, 2] x nominal, latent heat <= 3 x nominal, ambient +-10 K at arming, thermocouple lag <= 4 s, hard trip 130 C, T1-T2 disagreement 12 K (25 K while changing fast), stall: > A_f,max + 5 K for 8 s with < 15 % of the stroke moved | ASSUMED | tune on the bench (Gate 4 of BENCH_PROTOCOL.md) |
| Structural plant | 21 parallel elements, temperature spread 1-5 K (Gaussian) plus optional one-sided tail (25 % of elements, to 20 K) | ASSUMED | strain-temperature curve of the real wire |
| Other alloys | MEFF constants varied one at a time (h 100-1400, beta 4.5-8, eps_L 3-6 %, R 120-220, T0 +13/-10 K); thermal and electrical data kept | ASSUMED | datasheet of the actual alloy |

What the verification does NOT cover: the real wire (end effects at the clamps, non-uniform temperature, functional fatigue, creep), a different alloy's thermal data, the wing (coupled stiffness), the two-way effect, MATLAB/Simulink.
