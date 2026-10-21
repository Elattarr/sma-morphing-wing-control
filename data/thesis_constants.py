"""Wire and thermal constants copied from the thesis (literature values, not measured)."""
WIRE = dict(
    d_mm=2.0, L_mm=200.0, V_mm3=628.0, S_m2=1.26e-3,
    rho=6450.0, C=320.0, k=18.0,
    latent_J_per_mm3=0.156,
    rho_E_A=1.0e-6, rho_E_M=0.8e-6,      # ohm m
    h=81.9,                               # W/m^2/K (used); 5 also considered
    R_M=0.051, R_A=0.064,                 # ohm
)
T_AMBIENT = 22.0
AS_FREE, AF_FREE = 26.0, 52.0             # C, free wire (derived from MEFF constants)
FROZEN_STRAIN = 0.0327                    # sqrt(2/3) * 4 %
CONTRACTION_MM = FROZEN_STRAIN * WIRE["L_mm"]   # about 6.54 mm
HEATING_CASES = [(2.5, 35, 69, 173), (5.0, 26, 36, 182), (10.0, 19, 20, 199)]  # t s, I A, P W, E J
