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
# MEFF constants (ANSYS Shape Memory Effect) and free-wire transformation temperatures
# Auricchio & Petrini 2004: window Af - As = h*eps_L/beta (NOT h/beta)
MEFF = dict(h=145.0, T0=22.0, R=167.0, beta=5.6, eps_L=0.04)   # MPa, C, MPa, MPa/C, -
AF_FREE = MEFF["T0"] + MEFF["R"] / MEFF["beta"]                                    # 51.8 C
AS_FREE = MEFF["T0"] + (MEFF["R"] - MEFF["h"] * MEFF["eps_L"]) / MEFF["beta"]       # 50.8 C
FROZEN_STRAIN = 0.0327                    # sqrt(2/3) * 4 %
CONTRACTION_MM = FROZEN_STRAIN * WIRE["L_mm"]   # about 6.54 mm
HEATING_CASES = [(2.5, 35, 69, 173), (5.0, 26, 36, 182), (10.0, 19, 20, 199)]  # t s, I A, P W, E J
