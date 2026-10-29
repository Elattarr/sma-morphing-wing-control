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
