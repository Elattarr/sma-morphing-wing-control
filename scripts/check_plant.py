"""Checks of the standalone plant against ANSYS and the thesis heating cases.
Run from the repo root:  python scripts/check_plant.py"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "src"))
from sma_control import SMAParams, SMAWire

p = SMAParams()
m = p.meff
print(f"MEFF: As = {m.As:.2f} C, Af = {m.Af:.2f} C, window = {m.window:.2f} K")
print(f"wire: m = {p.mass*1e3:.2f} g, R_M = {p.resistance(0)*1e3:.1f} mohm, R_A = {p.resistance(1)*1e3:.1f} mohm, "
      f"m*c = {p.heat_capacity:.2f} J/K, latent = {p.latent_total:.1f} J, tau_cool = {p.cooling_tau:.1f} s")

# ---- 1. recovery map vs ANSYS 3D standalone wire (heating after force removal) ----
d = pd.read_csv(root / "data" / "ansys_3d_wire_merged.csv")
s = d[d.step.isin([4, 5]) & (d["T"] > 49)].sort_values("T")
r_ansys = 1.0 - s.p_max.values / s.p_max.values.max()      # recovered fraction
w = SMAWire(p)
# temperature error at equal recovery (more meaningful than amplitude error on a 1 K step)
inside = (r_ansys > 0.05) & (r_ansys < 0.95)
dT = np.array([w.map.T_of_r(r) for r in r_ansys[inside]]) - s["T"].values[inside]
print(f"\n[1] map vs ANSYS: temperature offset at equal recovery (5-95 %): mean {dT.mean():+.2f} K, max |{np.abs(dT).max():.2f}| K")
print(f"    ANSYS recovery starts ~{s['T'].values[r_ansys > 0.001].min():.2f} C, ends ~{s['T'].values[r_ansys >= 0.999].min():.2f} C"
      f"  (model: {p.meff.As:.2f} -> {p.meff.Af:.2f} C)")

# ---- 2. heating energy to 73 C vs thesis cases (R5) -------------------------------
print("\n[2] constant-current heating 22 -> 73 C (thesis: t, I, E = 2.5 s/35 A/173 J, 5 s/26 A/182 J, 10 s/19 A/199 J)")
for I, t_th, E_th in [(35, 2.5, 173), (26, 5.0, 182), (19, 10.0, 199)]:
    w = SMAWire(SMAParams(Imax=50))
    t = 0.0
    while w.T < 73.0 and t < 200:
        w.step(I, 0.001); t += 0.001
    print(f"    I = {I:2d} A: model t = {t:5.2f} s, E = {w.E_joule:6.1f} J   | thesis t = {t_th} s, E = {E_th} J")

# ---- 3. holding current at the end of the window ----------------------------------
for T in (m.As, m.Af):
    Ih = np.sqrt(p.h * p.surface * (T - p.T_inf) / p.resistance(1.0 if T >= m.Af else 0.5))
    print(f"\n[3] steady holding current at {T:.1f} C: ~{Ih:.1f} A, power {p.h*p.surface*(T-p.T_inf):.2f} W")
