"""Heating-only (up-stroke) controller on the thesis wing wire: 2 mm, 200 mm, structure only (no bias spring).

Plant: the wire pulls the wing structure, stiffness at the wire 6.4 N/mm, residual preload 7 N (2.2 MPa) after the ANSYS un-loading
(thesis Up 40 mm run, `deflection vs temp pos40` + `force vs temp position up 40mm`). MEFF is one-way: on cooling nothing happens, so
the wire cannot extend by itself; the controller only drives the wire upward (target r rises), holds without power, and a release
(power off + external reset) is outside the controller. No homing (it would use the stroke up). Windows come from the model
(A_s 51.2, A_f 54.1 C), the controller adapts its temperature limit when the real window is wider or higher.

TE deflection is mapped linearly from the wire recovery r:  TE = r * 42.2 mm  [ASSUMED: linear, from ANSYS force/deflection 1.15 N/mm].

usage:  python scripts/evaluate_heatonly.py steps|ramps|disturb|robust|sensors|struct [figure_dir]
"""
import os, sys
os.environ.update(SMA_D_MM="2", SMA_IMAX="35", SMA_SPRING_K="6.36", SMA_SIG="2.2", SMA_NO_HOMING="1")
sys.path.insert(0, os.path.dirname(__file__))
from dataclasses import replace
from pathlib import Path
import numpy as np
import evaluate_controller as ev
from sma_control.cycle_controller import CycleController, cfg_heat_only
from sma_control.spring import BiasSpring

TE_MM = 42.2                                   # thesis Up 40 mm, final TE deflection [mm] for r = 1
UP = (0.25, 0.5, 0.75, 1.0)
K_CTL = [6.36]                                 # structure stiffness the CONTROLLER assumes [N/mm] (the plant keeps 6.36)


HEAT_CFG = {}                              # K below A_s where the pre-heat stops (one-way: an overshoot cannot be undone)


def _make_ctl(sigma_ctl, cfg, p=None):
    p = p or ev.NOM
    return CycleController(p, sigma_ctl, cfg, K_spring=BiasSpring(sigma_ctl, K_CTL[0]).K(p))


ev.make_cfg = lambda: cfg_heat_only(ev.NOM, Ts=ev.TS, overshoot_exit=ev.OVERSHOOT_EXIT, **HEAT_CFG)
ev.make_ctl = _make_ctl                         # the harness calls make_ctl by module attribute


def sec_steps(fig=None):
    print("== HEAT-ONLY: STEP RESPONSES from r = 0 (2 mm wire, structure stiffness 6.4 N/mm, preload 7 N) ==")
    print(f"{'step':>12s} {'TE [mm]':>8s} {'delay':>6s} {'rise':>6s} {'settle':>7s} {'OS %':>6s} {'sse %':>7s} {'Tmax':>6s} {'Imax':>5s} {'E[J]':>6s}")
    runs = {}
    for r1 in UP:
        L = ev.simulate(lambda t, r1=r1: r1, 90.0, xi0=1.0, T0=22.0)
        m = ev.step_metrics(L, 0.0, r1)
        runs[r1] = L
        print(f"  0.00→{r1:.2f} {r1 * TE_MM:8.1f} {m['delay']:6.1f} {m['rise']:6.1f} {m['settle']:7.1f} {m['os']:6.2f} {m['sse']:+7.2f} {m['Tmax']:6.1f} {m['Imax']:5.1f} {m['E']:6.0f}")
    if fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        ev.sty(plt)
        f, ax = plt.subplots(1, 3, figsize=(15, 4.6))
        for r1, L in runs.items():
            ax[0].plot(L[:, 0], L[:, 2] * TE_MM, lw=2, label=f"target {r1 * TE_MM:.1f} mm"); ax[0].axhline(r1 * TE_MM, color=ev.GR, lw=1, ls="--")
            ax[1].plot(L[:, 0], L[:, 1], lw=1.6); ax[2].plot(L[:, 0], L[:, 4], lw=1.4)
        for a in ax:
            a.set_xlim(0, 30); a.set_xlabel("time [s]")
        ax[0].set_ylabel("TE deflection [mm]"); ax[1].set_ylabel("wire temperature [°C]"); ax[2].set_ylabel("current [A]"); ax[0].legend(frameon=False, fontsize=8.5)
        f.tight_layout(); f.savefig(Path(fig) / "heatonly_steps.png", dpi=150, facecolor=ev.BG)


def sec_ramps(fig=None):
    print("== HEAT-ONLY: RAMP TRACKING (reference 0 -> 1) ==")
    print(f"{'case':28s} {'RMS err %':>9s} {'RMS [mm]':>9s} {'max err %':>9s} {'best lag s':>10s} {'Tmax':>6s} {'E[J]':>6s}")
    for T in (20.0, 60.0, 180.0, 600.0):
        L = ev.simulate(lambda t, T=T: float(np.clip(t / T, 0, 1)), T + 40.0, xi0=1.0, T0=22.0)
        m = ev.tracking_metrics(L)
        print(f"up ramp 0→1 in {T:4.0f} s          {m['rms']:9.2f} {m['rms'] * TE_MM / 100:9.2f} {m['emax']:9.2f} {m['lag']:10.0f} {L[:, 1].max():6.1f} {L[-1, 6]:6.0f}")


def sec_disturb(fig=None):
    print("== HEAT-ONLY: DISTURBANCES while moving to r = 0.5 and to 1.0 ==")
    print(f"{'case':28s} {'target':>6s} {'time to target':>14s} {'final err %':>11s} {'Tmax':>6s} {'E[J]':>6s}")
    cases = [("nominal", {}), ("air flow x3 from 0 s", dict(h_scale=lambda t: 3.0)), ("air flow x0.5 from 0 s", dict(h_scale=lambda t: 0.5)),
             ("ambient -10 K", dict(T_amb=lambda t: 12.0, T0=12.0)), ("ambient +10 K", dict(T_amb=lambda t: 32.0, T0=32.0)),
             ("air flow x3 from 4 s", dict(h_scale=lambda t: 3.0 if t > 4.0 else 1.0))]
    for tgt in (0.5, 1.0):
        for name, kw in cases:
            L = ev.simulate(lambda t, tgt=tgt: tgt, 120.0, xi0=1.0, **{**dict(T0=22.0), **kw})
            tt = ev.time_to(L, tgt)
            print(f"{name:28s} {tgt:6.1f} {tt:14.1f} {100 * (L[L[:, 0] > 100, 2].mean() - tgt):+11.2f} {L[:, 1].max():6.1f} {L[-1, 6]:6.0f}")


def cases():
    lam = 0.0727
    m = ev.NOM.meff
    return [("nominal", {}), ("Lambda = 46 J", dict(p_true=replace(ev.NOM, latent_J_per_mm3=lam))), ("c = 320", dict(p_true=replace(ev.NOM, c=320.0))),
            ("h_c = 5", dict(p_true=replace(ev.NOM, h=5.0))), ("h_c = 82", dict(p_true=replace(ev.NOM, h=82.0))),
            ("rho reversed", dict(p_true=replace(ev.NOM, rho_eM=1.0e-6, rho_eA=0.87e-6))),
            ("rho +20 %", dict(p_true=replace(ev.NOM, rho_eM=0.96e-6, rho_eA=1.2e-6))), ("rho -20 %", dict(p_true=replace(ev.NOM, rho_eM=0.64e-6, rho_eA=0.8e-6))),
            ("preload ctl 2 -> 30 MPa", dict(sigma_ctl=30.0)), ("preload true 30 MPa", dict(sigma_true=30.0, sigma_ctl=2.2)),
            ("spread 2 K", dict(struct=dict(spread_K=2.0))), ("spread 5 K", dict(struct=dict(spread_K=5.0))),
            ("spread 2 K + 25 % tail to 20 K (ANSYS-like)", dict(struct=dict(spread_K=2.0, tail_frac=0.25, tail_K=20.0))),
            ("ambient 15 C", dict(T_amb=15.0)), ("ambient 32 C", dict(T_amb=32.0))]


def sec_robust(fig=None):
    print("== HEAT-ONLY: PARAMETER ERRORS, sequence r = 0.25 -> 0.5 -> 0.75 -> 1.0 (controller from nominal numbers) ==")
    print(f"{'case':44s} {'reached':>9s} {'max hold err %':>14s} {'[mm]':>6s} {'Tmax':>6s} {'time [s]':>9s} {'E [J]':>7s}")
    rows = []
    for name, kw in cases():
        L, seg, home = ev.sequence(targets=UP, hold_s=10.0, t_max=1500.0, label=name, **kw)
        s = ev.seq_summary(L, seg, home, n_targets=4)
        rows.append((name, s))
        print(f"{name:44s} {s['reached']:>5d} of 4 {s['max_err']:14.2f} {s['max_err'] * TE_MM / 100:6.2f} {s['Tmax']:6.1f} {s['t_total']:9.0f} {s['E']:7.0f}")
    for K in (0.5, 2.0):                                                    # the controller's stiffness is wrong by a factor
        K_CTL[0] = 6.36 * K
        L, seg, home = ev.sequence(targets=UP, hold_s=10.0, t_max=1500.0, label=f"K ctl x{K}")
        s = ev.seq_summary(L, seg, home, n_targets=4)
        print(f"{'controller stiffness x' + str(K):44s} {s['reached']:>5d} of 4 {s['max_err']:14.2f} {s['max_err'] * TE_MM / 100:6.2f} {s['Tmax']:6.1f} {s['t_total']:9.0f} {s['E']:7.0f}")
    K_CTL[0] = 6.36
    if fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        ev.sty(plt)
        L, seg, home = ev.sequence(targets=UP, hold_s=10.0, t_max=600.0, label="plot")
        f, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        ax[0].plot(L[:, 0], L[:, 3] * TE_MM, color=ev.OR, ls="--", label="target"); ax[0].plot(L[:, 0], L[:, 2] * TE_MM, color=ev.INK, lw=2, label="TE deflection")
        ax[1].plot(L[:, 0], L[:, 1], color=ev.BL); ax[2].plot(L[:, 0], L[:, 4], color=ev.BL)
        ax[0].set_ylabel("TE deflection [mm]"); ax[1].set_ylabel("T [°C]"); ax[2].set_ylabel("I [A]"); ax[2].set_xlabel("time [s]"); ax[0].legend(frameon=False)
        f.tight_layout(); f.savefig(Path(fig) / "heatonly_sequence.png", dpi=150, facecolor=ev.BG)


def sec_sensors(fig=None):
    print("== HEAT-ONLY: SENSORS, sequence r = 0.25 -> 0.5 -> 0.75 -> 1.0 ==")
    print(f"{'case':36s} {'reached':>9s} {'max hold err %':>14s} {'[mm]':>6s} {'Tmax':>6s} {'time [s]':>9s}")
    for name, kw in [("nominal (0.1 K, 0.02 mm)", {}), ("T noise 0.5 K", dict(Tn=0.5)), ("displacement noise 0.1 mm", dict(dn=0.1)),
                     ("displacement quantised 0.05 mm", dict(quant=0.05)), ("thermocouple lag 1 s", dict(T_lag=1.0)), ("thermocouple lag 3 s", dict(T_lag=3.0)),
                     ("all: 0.3 K, 0.05 mm, q 0.02, lag 1 s", dict(Tn=0.3, dn=0.05, quant=0.02, T_lag=1.0))]:
        L, seg, home = ev.sequence(targets=UP, hold_s=10.0, t_max=1500.0, sens=ev.Sensors(**kw), label=name)
        s = ev.seq_summary(L, seg, home, n_targets=4)
        print(f"{name:36s} {s['reached']:>5d} of 4 {s['max_err']:14.2f} {s['max_err'] * TE_MM / 100:6.2f} {s['Tmax']:6.1f} {s['t_total']:9.0f}")


T_TRIP = 85.0                                  # K  independent hard trip for the embedded wire (controller cap 75 C): [ASSUMED] set from the resin Tg


_new_guard = ev.new_guard


def _guard(cfg, sens=None, **kw):
    g = _new_guard(cfg, sens=sens, **kw)
    if g is not None:
        g.mon.c.T_trip = T_TRIP
    return g


ev.new_guard = _guard


def sec_safety(fig=None):
    print(f"== HEAT-ONLY: SAFETY, injected faults during the staircase (controller cap 75 C, monitor hard trip {T_TRIP:.0f} C) ==")
    cfg = ev.make_cfg()
    faults = [("thermocouple open (reads 0 C)", "tc_open"), ("thermocouple frozen", "tc_stuck"), ("thermocouple reads 30 K low", "tc_offset"),
              ("position sensor frozen", "d_stuck"), ("controller hangs, current frozen", "hang"), ("power stage stuck on", "stuck_on"),
              ("wire contact fails (R x3)", "wire_open")]
    print(f"{'fault':36s} {'at [s]':>6s} | {'no monitor: Tmax after':>22s} | {'stopped by':>22s} {'t [s]':>7s} {'Tmax after':>10s}  result")
    n = n_ok = 0
    for use_T2 in (True, False):
        for name, f in (faults if use_T2 else faults[:3]):
            for t_f in (3.0, 20.0, 40.0):
                res = []
                for protect in (False, True):
                    g = ev.Guard(ev.NOM, cfg.Ts, fault=f, t_fault=t_f, label=name, protect=protect, use_T2=use_T2)
                    g.mon.c.T_trip = T_TRIP
                    L, seg, home = ev.sequence(targets=UP, hold_s=10.0, t_max=300.0, guard=g, sens=ev.Sensors(T0=22.0))
                    by = g.mon.reason if (protect and g.mon.tripped) else None
                    res.append((g.Tmax_after, by, g.t))
                (Tn, _, _), (Tp, by, tt) = res
                ok = (by is not None and Tp < T_TRIP + 10.0) or Tn < T_TRIP
                n, n_ok = n + 1, n_ok + ok
                print(f"{name + ('' if use_T2 else ' [no T2]'):36s} {t_f:6.0f} | {Tn:22.1f} | {str(by):>22s} {tt:7.1f} {Tp:10.1f}  {'ok' if ok else 'NOT STOPPED'}")
    print(f"{n_ok} of {n} fault runs ended below the trip temperature + 10 K (or never exceeded the trip without the monitor); fail-safe = power off")


SECTIONS = dict(safety=sec_safety, steps=sec_steps, ramps=sec_ramps, disturb=sec_disturb, robust=sec_robust, sensors=sec_sensors)

if __name__ == "__main__":
    SECTIONS[sys.argv[1]](sys.argv[2] if len(sys.argv) > 2 else None)
    print("safety monitor: " + ("no trips in any run" if not ev.TRIPS else "TRIPS " + "; ".join(f"{a}: {b} at {c:.0f} s" for a, b, c in ev.TRIPS)))
