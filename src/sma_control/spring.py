"""One-way wire with a BIAS SPRING: the wire stress is no longer constant, it follows the displacement.

Mechanics (wire in series with a linear spring, spring stretched further when the wire contracts):

    sigma = sigma_ext + K (eps_ref - eps)               spring law, K = k L / A  [MPa per unit strain]
    eps   = sigma / E(xi) + eps_tr xi                   wire strain (E(xi) linear mixture, as in loading.py)

    eps_ref = sigma_ext / E_M + eps_tr                  strain of the fully extended wire (xi = 1), where sigma = sigma_ext

Eliminating eps gives sigma(xi) in closed form (SpringStress.sigma). K = 0 is the constant-stress wire of cycle_plant.py.

Phase law with the coupling: the ratchet limits depend on sigma, and sigma depends on xi, so each step solves

    xi_new = clip(xi_old, xi_f(sigma(xi_new), T), xi_r(sigma(xi_new), T))

g(x) = x - clip(...) is strictly increasing in x (sigma falls as xi rises, both limits fall with it), so the root is unique.

What this adds to the constant-stress wire (all follow from the two equations above, nothing is fitted):
  * heating contracts the wire, the spring force rises, so the recovery temperature rises during recovery:
        T_rec(xi) = T0 + (sqrt(2/3) sigma(xi) + R/(1-kappa m) - h eps_L xi) / beta
    the recovery window is wider than the free-wire 1 K by about sqrt(2/3) K Delta_eps / beta
  * cooling: the martensite is soft and the spring pulls it out, but only if the forward transformation can start:
        T_fwd(xi) = T0 + (sqrt(2/3) sigma(xi) - R/(1+kappa m) - h eps_L xi) / beta
    must stay above the ambient temperature all the way to xi = 1. At 22 C this needs sigma_ext >= about 211 MPa
    (sigma_ext is the SMALLEST stress the wire sees, because the spring force is lowest when the wire is extended).
"""
from dataclasses import dataclass
import math
from .params import SMAParams
from .loading import UniaxialMEFF
from .cycle_plant import SMAWireCycle


@dataclass(frozen=True)
class BiasSpring:
    sigma_ext: float = 250.0     # MPa  wire stress with the wire fully extended (spring preload / A)      [DESIGN]
    k: float = 20.0              # N/mm spring rate                                                         [DESIGN]

    def K(self, p: SMAParams) -> float:
        """MPa of wire stress per unit strain: k [N/m] * L [m] / A [m2] / 1e6"""
        return self.k * 1e3 * p.L / p.A / 1e6


class SpringStress:
    """closed-form sigma(xi), strain(xi), its inverse, and the two window temperatures T_rec(xi), T_fwd(xi)"""

    def __init__(self, u: UniaxialMEFF, sigma_ext: float, K: float):
        self.u, self.s0, self.K = u, sigma_ext, K
        self.eps_ref = sigma_ext / u.E_M + u.eps_tr
        self._c0 = sigma_ext + K * self.eps_ref
        M = u.meff
        self._q_rec = M.R / (1.0 - u.kappa * M.m)
        self._q_fwd = M.R / (1.0 + u.kappa * M.m)
        self._a = math.sqrt(2.0 / 3.0)

    def strain(self, xi: float) -> float:
        E = self.u.E(xi)
        return (self._c0 / E + self.u.eps_tr * xi) / (1.0 + self.K / E)

    def sigma(self, xi: float) -> float:
        return self.s0 + self.K * (self.eps_ref - self.strain(xi))

    def xi_from_strain(self, eps: float, n: int = 8) -> float:
        """inverse of strain(xi) (fixed point; contraction factor about K/E_A * eps_tr, i.e. about 1e-2)"""
        u, xi = self.u, 0.5
        for _ in range(n):
            xi = min(max((eps - self.sigma(xi) / u.E(xi)) / u.eps_tr, 0.0), 1.0)
        return xi

    def T_rec(self, xi: float) -> float:
        """temperature at which the reverse limit equals xi (recovery passes through xi)"""
        M = self.u.meff
        return M.T0 + (self._a * self.sigma(xi) + self._q_rec - M.h * M.eps_L * xi) / M.beta

    def T_fwd(self, xi: float) -> float:
        """temperature at which the forward limit equals xi (re-extension passes through xi)"""
        M = self.u.meff
        return M.T0 + (self._a * self.sigma(xi) - self._q_fwd - M.h * M.eps_L * xi) / M.beta


def _illinois(f, a, b, fa, fb, tol=1e-12, maxit=60):
    """bracketed regula falsi (Illinois variant) for a monotone f with f(a) and f(b) of opposite sign"""
    side, x = 0, 0.5 * (a + b)
    for _ in range(maxit):
        x = (a * fb - b * fa) / (fb - fa)
        fx = f(x)
        if abs(fx) < tol or abs(b - a) < 1e-14:
            return x
        if fx * fb > 0:
            b, fb = x, fx
            if side == -1:
                fa *= 0.5
            side = -1
        else:
            a, fa = x, fx
            if side == 1:
                fb *= 0.5
            side = 1
    return x


class SMAWireSpring(SMAWireCycle):
    """same energy balance and ratchet as SMAWireCycle, but sigma = sigma(xi) from the bias spring (no external stress input)"""

    def __init__(self, p: SMAParams | None = None, spring: BiasSpring | None = None, T0: float = 22.0, xi0: float = 1.0):
        super().__init__(p, T0, xi0)
        self.spring = spring or BiasSpring()
        self.st = SpringStress(self.u, self.spring.sigma_ext, self.spring.K(self.p))

    @property
    def sigma(self) -> float:
        return self.st.sigma(self.xi)

    # ---- phase law with the stress coupling (the sigma argument is ignored: the spring sets it)
    def xi_next(self, xi: float, sigma: float, T: float) -> float:
        u, st = self.u, self.st

        def c(x):
            s = st.sigma(x)
            return min(max(min(max(xi, u.xi_f(s, T)), u.xi_r(s, T)), 0.0), 1.0)

        if abs(c(xi) - xi) < 1e-13:                       # no transformation at this state
            return xi
        g0, g1 = -c(0.0), 1.0 - c(1.0)
        if g0 >= 0.0:
            return 0.0
        if g1 <= 0.0:
            return 1.0
        return _illinois(lambda x: x - c(x), 0.0, 1.0, g0, g1)

    def _solve_T(self, Tf: float, sigma: float) -> float:
        mc, mL = self.p.heat_capacity, self.p.latent_total
        lo, hi = Tf - mL / mc, Tf + mL / mc
        g = lambda T: mc * (Tf - T) + mL * (self.xi_next(self.xi, sigma, T) - self.xi)
        return _illinois(g, lo, hi, g(lo), g(hi), tol=1e-10)

    def step(self, I: float, dt: float, T_inf: float | None = None, h_scale: float = 1.0):
        return super().step(I, dt, 0.0, T_inf, h_scale)

    def disp_mm(self, sigma: float = 0.0) -> float:
        return self.p.L * 1e3 * self.st.strain(self.xi)
