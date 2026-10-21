"""Parameters of the standalone (free, one-way) SMA wire plant.

Every value carries a tag:
  [THESIS]   taken from the thesis (MEFF constants, ANSYS results, literature values quoted there)
  [ANSYS]    read from the ANSYS standalone-wire results in data/
  [DESIGN]   a design choice of the controller/hardware, not a material property
  [ASSUMED]  placeholder that bench data must replace (see README, list of assumptions)
"""
from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class MEFF:
    """ANSYS Shape Memory Effect (Auricchio-Petrini) constants [THESIS]."""
    h: float = 145.0        # hardening parameter [MPa]
    T0: float = 22.0        # phase-balance temperature [C]
    R: float = 167.0        # elastic-domain radius (deviatoric) [MPa]
    beta: float = 5.6       # temperature scaling [MPa/C]
    eps_L: float = 0.04     # max transformation strain (equivalent)

    @property
    def Af(self) -> float:                      # end of recovery, zero stress
        return self.T0 + self.R / self.beta

    @property
    def As(self) -> float:                      # start of recovery, zero stress
        return self.T0 + (self.R - self.h * self.eps_L) / self.beta

    @property
    def window(self) -> float:                  # Af - As = h*eps_L/beta
        return self.h * self.eps_L / self.beta

    @property
    def stress_shift(self) -> float:            # uniaxial stress shift of the transformation temperatures [MPa/K]
        return math.sqrt(1.5) * self.beta


@dataclass
class SMAParams:
    # geometry [THESIS]
    d: float = 2.0e-3                 # diameter [m]
    L: float = 0.20                   # active length [m]
    # material / thermal [THESIS, literature values quoted in the thesis]
    rho: float = 6450.0               # density [kg/m3]
    c: float = 320.0                  # specific heat [J/kg/K]
    latent_J_per_mm3: float = 0.156   # latent heat of M->A [J/mm3]  (= 24.2 kJ/kg)
    rho_eM: float = 0.8e-6            # resistivity martensite [ohm m]
    rho_eA: float = 1.0e-6            # resistivity austenite [ohm m]
    h: float = 81.9                   # convection coefficient [W/m2/K] (thesis: 81.9 used, 5 considered)
    T_inf: float = 22.0               # ambient [C]
    # recovery map
    meff: MEFF = MEFF()
    full_stroke_mm: float = 6.64      # [ANSYS] contraction at zero stress (3.32 % of 200 mm)
    tail_fraction: float = 0.0        # [ASSUMED] share of the stroke recovered in a slow tail above Af
    tail_width: float = 40.0          # [ASSUMED] width of that tail [K]
    T_shift: float = 0.0              # [ASSUMED] shift of the whole map [K] (ANSYS runs ~0.1 K above the formula)
    # hardware / sensing
    Imax: float = 35.0                # [THESIS R6] peak current available [A]
    noise_mm: float = 0.02            # [ASSUMED] position-sensor noise, 1 sigma [mm]
    quant_mm: float = 0.0             # [ASSUMED] sensor resolution [mm], 0 = off

    # ---- derived ---------------------------------------------------------
    @property
    def A(self) -> float:
        return math.pi * self.d ** 2 / 4

    @property
    def V(self) -> float:
        return self.A * self.L

    @property
    def mass(self) -> float:
        return self.rho * self.V

    @property
    def surface(self) -> float:
        return math.pi * self.d * self.L

    @property
    def latent_J_per_kg(self) -> float:
        return self.latent_J_per_mm3 * 1e9 / self.rho

    @property
    def heat_capacity(self) -> float:           # m*c [J/K]
        return self.mass * self.c

    @property
    def latent_total(self) -> float:            # m*Lambda [J]
        return self.mass * self.latent_J_per_kg

    @property
    def cooling_tau(self) -> float:             # m*c/(h*S) [s]
        return self.heat_capacity / (self.h * self.surface)

    def resistance(self, r: float) -> float:
        """wire resistance at recovered fraction r (0 martensite, 1 austenite) [ohm]"""
        rho_e = self.rho_eM + r * (self.rho_eA - self.rho_eM)
        return rho_e * self.L / self.A


class RecoveryMap:
    """Recovered fraction r (0 = frozen martensite, 1 = fully recovered) vs temperature
    at zero stress: piecewise linear.  MEFF is linear between As and Af (linear hardening);
    an optional slow tail above Af mimics real wires (assumed, to be fitted to bench data)."""

    def __init__(self, p: SMAParams):
        As, Af = p.meff.As + p.T_shift, p.meff.Af + p.T_shift
        if p.tail_fraction > 0:
            self.T = np.array([As, Af, Af + p.tail_width])
            self.r = np.array([0.0, 1.0 - p.tail_fraction, 1.0])
        else:
            self.T = np.array([As, Af])
            self.r = np.array([0.0, 1.0])

    def r_of_T(self, T: float) -> float:
        return float(np.interp(T, self.T, self.r, left=0.0, right=1.0))

    def T_of_r(self, r: float) -> float:
        return float(np.interp(r, self.r, self.T))
