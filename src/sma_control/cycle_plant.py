"""Full-cycle wire: heating AND cooling, with an external stress input and an optional two-way extension.

State: T (wire temperature), xi (martensite fraction, 0..1).   r = 1 - xi is the recovered fraction of plant.py.
Input: I (current), sigma_ext (axial stress from the load or bias spring, MPa), T_inf.

Phase law (same ratchet as loading.py):   xi_new = clip(xi, lo(sigma, T), hi(sigma, T))
    hi = xi_r(sigma, T)   reverse limit  -> recovery on heating
    lo = xi_f(sigma, T)   forward limit  -> re-extension on cooling, only when sigma is large enough
Energy (latent heat absorbed on recovery, released on forward transformation):
    m c dT/dt = I^2 R(r) - h S (T - T_inf) + m*Lambda * dxi/dt
solved implicitly each step by bisection (T and xi are coupled through the latent heat).

ONE-WAY MEFF only: at sigma = 0 the forward limit is far below the temperatures reached (about -8 C), so cooling changes
nothing. A stress sigma > 0 moves the forward start up to  T_start(sigma) = T0 + (sqrt(2/3) sigma - R/(1+kappa m))/beta.
The two-way extension (an assumption, not MEFF) lives separately in two_way.py.
"""
import numpy as np
from .params import SMAParams
from .loading import UniaxialMEFF


class SMAWireCycle:
    def __init__(self, p: SMAParams | None = None, T0: float = 22.0, xi0: float = 1.0):
        self.p = p or SMAParams()
        self.u = UniaxialMEFF(self.p.meff)
        self.xi = xi0
        self.T = T0
        self.E_joule = 0.0

    # ------------------------------------------------------------------ phase law
    def limits(self, sigma: float, T: float):
        return self.u.xi_f(sigma, T), self.u.xi_r(sigma, T)

    def xi_next(self, xi: float, sigma: float, T: float) -> float:
        lo, hi = self.limits(sigma, T)
        return min(max(min(max(xi, lo), hi), 0.0), 1.0)

    # ------------------------------------------------------------------ dynamics
    def step(self, I: float, dt: float, sigma: float = 0.0, T_inf: float | None = None, h_scale: float = 1.0):
        p = self.p
        Tinf = p.T_inf if T_inf is None else T_inf
        I = float(np.clip(I, 0.0, p.Imax))
        P = I * I * p.resistance(1.0 - self.xi)
        Q = h_scale * p.h * p.surface * (self.T - Tinf)
        mc, mL = p.heat_capacity, p.latent_total
        Tf = self.T + dt * (P - Q) / mc                       # temperature without latent heat
        xi_try = self.xi_next(self.xi, sigma, Tf)
        if abs(xi_try - self.xi) > 1e-12:                      # transformation active: latent heat couples T and xi
            T_new = self._solve_T(Tf, sigma)
            xi_new = self.xi_next(self.xi, sigma, T_new)
        else:
            T_new, xi_new = Tf, self.xi
        self.E_joule += P * dt
        self.T, self.xi = T_new, xi_new
        return self.T, self.xi

    def _solve_T(self, Tf: float, sigma: float) -> float:
        """temperature after the step when latent heat is exchanged: root of  mc (Tf - T) + mL (xi(T) - xi_old) = 0  (bisection)"""
        mc, mL = self.p.heat_capacity, self.p.latent_total
        lo, hi = Tf - mL / mc, Tf + mL / mc
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            g = mc * (Tf - mid) + mL * (self.xi_next(self.xi, sigma, mid) - self.xi)
            lo, hi = (mid, hi) if g > 0 else (lo, mid)
        return 0.5 * (lo + hi)

    def step_isothermal(self, sigma: float, T: float):
        """quasi-static loading with a perfect heat sink (the ANSYS loading steps are isothermal): T is imposed, no latent-heat coupling"""
        self.T = T
        self.xi = self.xi_next(self.xi, sigma, T)
        return self.T, self.xi

    # ------------------------------------------------------------------ outputs
    def disp_mm(self, sigma: float = 0.0) -> float:
        """end displacement, same convention as ANSYS defl_max (L * strain)"""
        return self.p.L * 1e3 * self.u.strain(sigma, self.xi)

    def T_forward_start(self, sigma: float) -> float:
        """temperature at which the forward limit reaches xi = 0 under external stress sigma"""
        M = self.p.meff
        return M.T0 + (np.sqrt(2.0 / 3.0) * sigma - M.R / (1.0 + self.u.kappa * M.m)) / M.beta
