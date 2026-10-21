"""Closed-loop step tests: python scripts/step_test.py"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sma_control import SMAParams, SMAWire
from sma_control.controller import Controller, CtrlCfg


def run(r_target, T_end=30.0, h_scale=1.0, T_amb=22.0, p=None, cfg=None, Tn=0.1, seed=1):
    p = p or SMAParams()
    w = SMAWire(p, T0=T_amb, seed=seed)
    cfg = cfg or CtrlCfg()
    ctl = Controller(p, cfg)
    rng = np.random.default_rng(seed + 100)
    n = int(T_end / cfg.Ts)
    log = np.zeros((n, 5))
    for k in range(n):
        T_meas = w.T + Tn * rng.standard_normal()
        I = ctl.step(r_target, T_meas, w.measure())
        w.step(I, cfg.Ts, T_inf=T_amb, h_scale=h_scale)
        log[k] = (k * cfg.Ts, w.T, w.r, I, ctl.Tcmd)
    return log, w


def metrics(log, r_target, w):
    t, T, r, I = log[:, 0], log[:, 1], log[:, 2], log[:, 3]
    ok = np.where(np.abs(r - r_target) <= 0.02)[0]
    ts = t[ok[0]] if len(ok) else np.nan
    return dict(t_settle=ts, r_final=r[-1], os_pct=100 * max(r.max() - r_target, 0) / max(r_target, 1e-9),
                Tmax=T.max(), Imax=I.max(), E=w.E_joule)


if __name__ == "__main__":
    for rt in (0.5, 1.0):
        for lab, kw in [("nominal", {}), ("h x0.2", dict(h_scale=0.2)), ("h x1 amb 35C", dict(T_amb=35.0)),
                        ("amb 15C", dict(T_amb=15.0))]:
            log, w = run(rt, **kw)
            m = metrics(log, rt, w)
            print(f"r_ref={rt:.1f} {lab:14s} ts={m['t_settle']:5.2f}s r_end={m['r_final']:.3f} OS={m['os_pct']:4.1f}% "
                  f"Tmax={m['Tmax']:.2f}C Imax={m['Imax']:.1f}A E={m['E']:.0f}J")
