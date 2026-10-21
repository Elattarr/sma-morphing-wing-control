# Plant assumptions and the measurement that replaces each

| Parameter | Value | Tag | Replace with |
|---|---|---|---|
| tail_fraction, tail_width | 0, 40 K | ASSUMED | strain-vs-temperature curve of the real wire at ~0 stress (heating, fixed small load): fit steep edge + tail |
| T_shift | 0 K | ASSUMED | same curve: measured As/Af vs MEFF 50.79/51.82 C |
| noise_mm, quant_mm | 0.02, 0 mm | ASSUMED | position-sensor datasheet / measured noise at rest |
| h (convection) | 81.9 W/m2K | THESIS | T-vs-time at fixed currents on the bench wire (h depends on airflow) |
| c, latent heat, rho_e(M/A) | literature | THESIS | DSC (c, latent); 4-wire resistance vs temperature |
| Imax | 35 A | THESIS R6 | actual supply/driver limit |
| MEFF h, R, beta, eps_L, T0 | 145, 167, 5.6, 0.04, 22 | THESIS | two strain-temperature runs at different stresses (beta), DSC |

Unresolved: ANSYS forward-transformation stress rise (2.6 vs 7.1 MPa predicted) - not used by this plant.
