"""Structurally different plant: a wire made of N parallel elements with a spread of transformation temperatures.

Why: the controller's model (cycle_controller.py) is the sharp single-element MEFF wire: one window of h eps_L / beta (1 K, or a few K
under a spring). Real wires have broader windows with long tails and hysteresis sub-loops. Testing the controller on its own model
only tests parameter errors; this plant tests a model error.

Model (a Preisach-type construction, nothing fitted):
  element i has its own state xi_i and a temperature offset d_i: it transforms as the MEFF element would at T - d_i
      xi_i,new = clip(xi_i,old, xi_f(sigma, T - d_i), xi_r(sigma, T - d_i))
  all elements carry the same stress, set by the bias spring from the MEAN state  xi = mean(xi_i)   (strain is shared)
      sigma = sigma(xi)   (spring.SpringStress)
  the wire state seen by everything else (resistance, latent heat, displacement) is the mean xi.
Each element keeps its own ratchet, so partial heating and cooling cycles give nested hysteresis sub-loops with return-point memory.

Offsets: a Gaussian part (standard deviation spread_K, placed at the quantiles of N elements, mean zero) plus an optional one-sided
tail: tail_frac of the elements with offsets spread evenly over 0..tail_K (long tail above A_f on heating, above the forward start on cooling).
Window width: about  h eps_L / beta + 2 * 1.64 * spread_K  between 5 % and 95 % transformed.
spread_K = 0 and tail_frac = 0 gives back SMAWireSpring exactly.
"""
from statistics import NormalDist
import numpy as np
from .params import SMAParams
from .spring import SMAWireSpring, BiasSpring, _illinois


def element_offsets(n: int, spread_K: float, tail_frac: float, tail_K: float) -> np.ndarray:
    n_tail = int(round(tail_frac * n))
    n_g = n - n_tail
    nd = NormalDist()
    g = np.array([spread_K * nd.inv_cdf((i + 0.5) / n_g) for i in range(n_g)]) if n_g else np.array([])
    t = np.linspace(0.0, tail_K, n_tail + 2)[1:-1] if n_tail else np.array([])
    return np.concatenate([g, t])


class SMAWireElements(SMAWireSpring):
    def __init__(self, p: SMAParams | None = None, spring: BiasSpring | None = None, T0: float = 22.0, xi0: float = 1.0,
                 n: int = 21, spread_K: float = 3.0, tail_frac: float = 0.0, tail_K: float = 0.0):
        super().__init__(p, spring, T0, xi0)
        self.offsets = element_offsets(n, spread_K, tail_frac, tail_K)
        self.xe = np.full(len(self.offsets), float(xi0))

    def _c(self, x: float, T: float, xe: np.ndarray) -> np.ndarray:
        s, u = self.st.sigma(x), self.u
        Te = T - self.offsets
        return np.clip(np.minimum(np.maximum(xe, u.xi_f(s, Te)), u.xi_r(s, Te)), 0.0, 1.0)

    def xi_next(self, xi: float, sigma: float, T: float) -> float:
        xe = self.xe
        f = lambda x: x - float(self._c(x, T, xe).mean())
        if abs(f(xi)) < 1e-13:                              # no transformation
            return xi
        g0, g1 = f(0.0), f(1.0)
        if g0 >= 0.0:
            return 0.0
        if g1 <= 0.0:
            return 1.0
        return _illinois(f, 0.0, 1.0, g0, g1)

    def step(self, I: float, dt: float, T_inf: float | None = None, h_scale: float = 1.0):
        xe_old = self.xe.copy()
        out = super().step(I, dt, T_inf, h_scale)
        self.xe = self._c(self.xi, self.T, xe_old)
        self.xi = float(self.xe.mean())
        return out
