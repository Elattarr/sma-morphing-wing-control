"""Standalone one-way SMA wire: current -> temperature -> recovery -> contraction.

States   : T (wire temperature), r (recovered fraction, monotonic)
Equations:
    m c dT/dt = I^2 R(r) - h S (T - T_inf) - m*Lambda * dr/dt        (energy balance, latent heat absorbed)
    r         = max(r_prev, r_map(T))                                 (one-way: no re-strain at zero stress)
    y         = full_stroke * r                                       (contraction [mm])
Cooling does not change r: with zero stress the MEFF model has no forward transformation.
"""
from dataclasses import dataclass
import numpy as np
from .params import SMAParams, RecoveryMap


@dataclass
class State:
    T: float
    r: float


class SMAWire:
    def __init__(self, p: SMAParams | None = None, T0: float | None = None, seed: int = 0):
        self.p = p or SMAParams()
        self.map = RecoveryMap(self.p)
        self.T = self.p.T_inf if T0 is None else T0
        self.r = 0.0
        self.rng = np.random.default_rng(seed)
        self.E_joule = 0.0          # electrical energy dissipated in the wire [J]

    # --------------------------------------------------------------------
    def step(self, I: float, dt: float, T_inf: float | None = None, h_scale: float = 1.0):
        """advance dt [s] with wire current I [A]; T_inf and h_scale are disturbances"""
        p = self.p
        Tinf = p.T_inf if T_inf is None else T_inf
        I = float(np.clip(I, 0.0, p.Imax))
        P = I * I * p.resistance(self.r)
        Q = h_scale * p.h * p.surface * (self.T - Tinf)
        mc = p.heat_capacity
        T_free = self.T + dt * (P - Q) / mc            # temperature without latent heat
        r_new = self.r
        T_new = T_free
        if self.map.r_of_T(T_free) > self.r:           # transformation proceeds -> absorbs latent heat
            lam = p.latent_total
            lo = T_free - lam * (1.0 - self.r) / mc    # lower bound: all remaining latent heat absorbed
            hi = T_free
            for _ in range(50):                        # root of  mc*(T_free - T) = lam*(r_map(T) - r)
                mid = 0.5 * (lo + hi)
                g = mc * (T_free - mid) - lam * (max(self.map.r_of_T(mid), self.r) - self.r)
                lo, hi = (mid, hi) if g > 0 else (lo, mid)
            T_new = 0.5 * (lo + hi)
            r_new = max(self.r, self.map.r_of_T(T_new))
        self.E_joule += P * dt
        self.T, self.r = T_new, r_new
        return self.T, self.r, self.stroke_mm()

    def stroke_mm(self) -> float:
        return self.p.full_stroke_mm * self.r

    def measure(self) -> float:
        p = self.p
        y = self.stroke_mm() + p.noise_mm * self.rng.standard_normal()
        if p.quant_mm > 0:
            y = round(y / p.quant_mm) * p.quant_mm
        return y

    def state(self) -> State:
        return State(self.T, self.r)
