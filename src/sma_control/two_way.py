"""Two-way extension of the cycle plant. NOT part of MEFF: a phenomenological add-on with two ASSUMED parameters.

A trained wire carries an internal stress sigma_int that adds to the external stress in the FORWARD limit only
(the reverse/recovery limit is unchanged), and the transformation it can drive is capped at xi_cap:

    lo = max( xi_f(sigma, T),  min( xi_f(sigma + sigma_int, T), xi_cap ) )

sigma_int follows from the chosen zero-stress forward start temperature Ms (T_start(sigma_int) = Ms).
Ms and xi_cap are not identified from data; they reproduce the qualitative effect only.
"""
from dataclasses import dataclass
from .params import SMAParams
from .cycle_plant import SMAWireCycle


@dataclass
class TwoWayParams:
    Ms: float = 35.0        # C, forward start at zero external stress (assumption)
    xi_cap: float = 0.3     # largest martensite fraction the internal stress alone creates (assumption)


class SMAWireTwoWay(SMAWireCycle):
    def __init__(self, p: SMAParams | None = None, T0: float = 22.0, xi0: float = 1.0, tw: TwoWayParams | None = None):
        super().__init__(p, T0, xi0)
        self.tw = tw or TwoWayParams()
        self.sig_int = self.u.sigma_start(self.tw.Ms)

    def limits(self, sigma: float, T: float):
        lo, hi = super().limits(sigma, T)
        lo = max(lo, min(self.u.xi_f(sigma + self.sig_int, T), self.tw.xi_cap))
        return lo, hi
