"""Map between the wire recovery r (0 = extended ... 1 = fully recovered) and the trailing-edge deflection of the wing in mm.

The controller works in r. The wing setpoint is given in mm, so a map is needed. Default: linear, 0 -> 0 mm, 1 -> -42.2 mm (thesis ANSYS,
Up 40 mm run; the sign is the ANSYS directional deformation, negative = TE down). [ASSUMED] linear: the thesis run gave a straight
force / deflection line of about 1.15 N/mm at the TE, but the deflection / recovery relation of the real wing (contact, CFRP, friction) is not
measured. Replace the points with a calibration (several (r, TE) pairs from the bench) with `TEMap(points)`; piecewise linear in between.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class TEMap:
    points: tuple = ((0.0, 0.0), (1.0, -42.2))       # (r, TE mm), r strictly increasing

    def __post_init__(self):
        r = [p[0] for p in self.points]
        if len(r) < 2 or any(b <= a for a, b in zip(r, r[1:])):
            raise ValueError("points need at least two entries with strictly increasing r")
        if abs(r[0]) > 1e-9 or abs(r[-1] - 1.0) > 1e-9:
            raise ValueError("the map must start at r = 0 and end at r = 1")

    def r_to_te(self, r):
        r_, te_ = zip(*self.points)
        return float(np.interp(np.clip(r, 0.0, 1.0), r_, te_))

    def te_to_r(self, te_mm):
        r_, te_ = zip(*self.points)
        te_ = np.array(te_)
        if te_[-1] < te_[0]:                        # decreasing in mm (TE down is negative): np.interp needs ascending x
            return float(np.interp(te_mm, te_[::-1], np.array(r_)[::-1]))
        return float(np.interp(te_mm, te_, r_))

    @property
    def te_range(self):
        te = [p[1] for p in self.points]
        return min(te), max(te)
