"""Uniaxial (1D) MEFF model with stress input: loading, unloading, heating at constant or zero stress.

State: xi in [0, 1]  (stored transformation strain / (sqrt(2/3)*eps_L)).
Input: axial stress sigma >= 0 [MPa] (tension) and temperature T [C].

Uniaxial reduction of the Auricchio-Petrini limit function F = |q| + kappa*m*q - R with
    q = sqrt(2/3)*sigma - [ beta*(T - T0) + h*eps_L*xi ]          (signed deviatoric distance, MPa)
    kappa = sqrt(2/27)
forward  (q > 0, austenite -> stress-induced martensite):  q = R / (1 + kappa*m)
reverse  (q < 0, martensite -> austenite)               :  q = -R / (1 - kappa*m)
which give the two ratcheting limits xi_f and xi_r below; the new value is clip(xi_prev, xi_f, xi_r).

At sigma = 0 this reduces to the recovery map used by plant.py (A_s, A_f, window) - see tests/test_loading.py.
"""
from dataclasses import dataclass
import math
from .params import MEFF


@dataclass
class UniaxialMEFF:
    meff: MEFF = MEFF()
    E_A: float = 51.7e3          # MPa, austenite modulus (ANSYS elastic slope at start of loading: 51.5 GPa)
    E_M: float = 46.0e3          # MPa, martensite modulus (ANSYS elastic slope after transformation: about 45.7 GPa)
    xi: float = 0.0

    # ---- constants ----------------------------------------------------------------
    @property
    def kappa(self) -> float:
        return math.sqrt(2.0 / 27.0)

    @property
    def eps_tr(self) -> float:                    # uniaxial transformation strain sqrt(2/3)*eps_L
        return math.sqrt(2.0 / 3.0) * self.meff.eps_L

    # ---- limits -------------------------------------------------------------------
    def xi_f(self, sigma: float, T: float) -> float:
        """forward (loading) limit: xi cannot be below this"""
        M = self.meff
        q0 = M.R / (1.0 + self.kappa * M.m)
        return (math.sqrt(2.0 / 3.0) * sigma - M.beta * (T - M.T0) - q0) / (M.h * M.eps_L)

    def xi_r(self, sigma: float, T: float) -> float:
        """reverse (recovery) limit: xi cannot be above this"""
        M = self.meff
        q0 = M.R / (1.0 - self.kappa * M.m)
        return (math.sqrt(2.0 / 3.0) * sigma - M.beta * (T - M.T0) + q0) / (M.h * M.eps_L)

    def step(self, sigma: float, T: float) -> float:
        xi = min(max(self.xi, self.xi_f(sigma, T)), self.xi_r(sigma, T))
        self.xi = min(max(xi, 0.0), 1.0)
        return self.xi

    # ---- mechanics ----------------------------------------------------------------
    def E(self, xi: float | None = None) -> float:
        xi = self.xi if xi is None else xi
        return self.E_A + xi * (self.E_M - self.E_A)

    def strain(self, sigma: float, xi: float | None = None) -> float:
        xi = self.xi if xi is None else xi
        return sigma / self.E(xi) + self.eps_tr * xi

    # ---- transformation temperatures at constant stress -------------------------------
    def A_s(self, sigma: float = 0.0) -> float:
        """temperature at which xi_r = 1 (recovery starts from a fully transformed wire)"""
        M = self.meff
        q0 = M.R / (1.0 - self.kappa * M.m)
        return M.T0 + (math.sqrt(2.0 / 3.0) * sigma + q0 - M.h * M.eps_L) / M.beta

    def A_f(self, sigma: float = 0.0) -> float:
        M = self.meff
        q0 = M.R / (1.0 - self.kappa * M.m)
        return M.T0 + (math.sqrt(2.0 / 3.0) * sigma + q0) / M.beta

    @property
    def shift_per_MPa(self) -> float:             # dA/dsigma [K/MPa]
        return math.sqrt(2.0 / 3.0) / self.meff.beta

    def sigma_start(self, T: float) -> float:     # forward onset stress at temperature T
        M = self.meff
        return (M.R / (1.0 + self.kappa * M.m) + M.beta * (T - M.T0)) / math.sqrt(2.0 / 3.0)

    def sigma_finish(self, T: float) -> float:
        return self.sigma_start(T) + self.meff.h * self.meff.eps_L / math.sqrt(2.0 / 3.0)
