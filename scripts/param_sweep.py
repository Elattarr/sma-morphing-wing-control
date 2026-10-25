"""Robustness sweep: true plant differs from the controller's nominal model.
python scripts/param_sweep.py   (controller always built from nominal SMAParams())"""
import sys
from dataclasses import replace
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sma_control import SMAParams, SMAWire
from sma_control.params import MEFF
from sma_control.controller import Controller, CtrlCfg

NOM = SMAParams()
LAM_CC = 0.0727            # J/mm3, Clausius-Clapeyron value (46 J)


def wide_window(width_K=12.0):
    """same A_f as nominal, window widened by the hardening parameter h"""
    m = NOM.meff
    h = width_K * m.beta / m.eps_L
    return replace(m, h=h, R=m.R)


def run(p_true, r_ref, T_amb=22.0, T_end=40.0, Tn=0.1, seed=1):
    cfg = CtrlCfg()
    w = SMAWire(p_true, T0=T_amb, seed=seed)
    ctl = Controller(NOM, cfg)
    rng = np.random.default_rng(seed + 100)
    n = int(T_end / cfg.Ts)
    t = np.arange(n) * cfg.Ts
    r = np.zeros(n); T = np.zeros(n); I = np.zeros(n)
    for k in range(n):
        Im = ctl.step(r_ref, w.T + Tn * rng.standard_normal(), w.measure())
        w.step(Im, cfg.Ts, T_inf=T_amb)
        r[k], T[k], I[k] = w.r, w.T, Im
    ok = np.where(np.abs(r - r_ref) <= 0.02)[0]
    return dict(ts=t[ok[0]] if len(ok) else np.nan, err=r[-1] - r_ref, os=100 * max(r.max() - r_ref, 0) / r_ref,
                Tmax=T.max(), Imax=I.max(), E=w.E_joule)


CASES = [
    ("nominal", NOM, 22.0),
    ("c = 320", replace(NOM, c=320.0), 22.0),
    ("Lambda = 46 J (CC)", replace(NOM, latent_J_per_mm3=LAM_CC), 22.0),
    ("h_c = 5", replace(NOM, h=5.0), 22.0),
    ("h_c = 20", replace(NOM, h=20.0), 22.0),
    ("h_c = 82", replace(NOM, h=82.0), 22.0),
    ("rho_e reversed (M 1.0 / A 0.87)", replace(NOM, rho_eM=1.0e-6, rho_eA=0.87e-6), 22.0),
    ("rho_e +20 %", replace(NOM, rho_eM=0.96e-6, rho_eA=1.2e-6), 22.0),
    ("rho_e -20 %", replace(NOM, rho_eM=0.64e-6, rho_eA=0.8e-6), 22.0),
    ("ambient 15 C", NOM, 15.0),
    ("ambient 35 C", NOM, 35.0),
    ("window 12 K (same A_f)", replace(NOM, meff=wide_window(12.0)), 22.0),
    ("worst mix: c 320, h 82, rho rev", replace(NOM, c=320.0, h=82.0, rho_eM=1.0e-6, rho_eA=0.87e-6), 22.0),
]

if __name__ == "__main__":
    print(f"{'case':36s} {'r_ref':>5s} {'ts[s]':>6s} {'err':>7s} {'OS%':>5s} {'Tmax':>6s} {'Imax':>5s} {'E[J]':>6s}")
    for name, p, Ta in CASES:
        for rr in (0.5, 1.0):
            m = run(p, rr, T_amb=Ta)
            print(f"{name:36s} {rr:5.1f} {m['ts']:6.2f} {m['err']:+7.3f} {m['os']:5.1f} {m['Tmax']:6.2f} {m['Imax']:5.1f} {m['E']:6.0f}")
