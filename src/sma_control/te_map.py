"""Map between the wire recovery r (0 = extended ... 1 = fully recovered, measured at the wire end) and the trailing-edge deflection in mm.

Default points: ANSYS Up 40 mm run (2026-10-30, heating 22 -> 80 C, docs/ansys_up40_recovery.csv). TE deflection is the ANSYS directional deformation
(negative = TE down) in the coordinates of that run: r = 0 is the PRELOADED state at 22 C (-6.63 mm, not 0), r = 1 is -42.24 mm, so the useful travel is 35.6 mm.
r is the recovery of the wire-end displacement (mean of min / max of the end-of-wire output, stroke 6.48 mm). The relation is concave: TE lags the wire end
by up to 13 % of the stroke at r = 0.2-0.4 (ANSYS data), so the old linear map was wrong in the middle of the stroke.
[ASSUMED] the same relation holds for another load case or wing; replace the points with a bench calibration, `TEMap(points)`, piecewise linear in between.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class TEMap:
    points: tuple = ((0.0, -6.632), (0.05, -7.155), (0.1, -8.144), (0.2, -10.537), (0.3, -13.252), (0.4, -16.449), (0.5, -20.748),
                      (0.6, -25.053), (0.7, -29.379), (0.8, -33.676), (0.9, -38.065), (0.95, -40.248), (1.0, -42.239))       # (r, TE mm), r strictly increasing

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
