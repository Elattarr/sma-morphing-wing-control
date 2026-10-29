"""Full one-way cycle of the standalone wire against ANSYS (steps 1-6).

ANSYS steps (data/ansys_3d_wire_merged.csv):
  1-3  load 0 -> 857 N at 22 C          4  unload while warming 22 -> 35 C
  5    heat 35 -> 52 C (recovery)       6  cool 52 -> 22 C
Our chain:
  loading / unloading : UniaxialMEFF (sigma, T) -> xi              (loading.py)
  hand-off            : r = 1 - xi/xi_stored,  stored displacement = L*eps_tr*xi
  heating             : plant (current -> T -> r) driven by the controller, start 35 C, T_amb 22 C
  cooling             : I = 0, plant cools with tau = m c /(h S)
  reload              : UniaxialMEFF again at 22 C  -> cycle 2
Displacement convention: d = y_stored - contraction  (same sign as ANSYS defl_max).

python scripts/full_cycle.py [figure_dir]
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sma_control import SMAParams, SMAWire
from sma_control.controller import Controller, CtrlCfg
from sma_control.loading import UniaxialMEFF
from sma_control.params import RecoveryMap

ROOT = Path(__file__).resolve().parents[1]
P0 = SMAParams()
A_MM2, L_MM = P0.A * 1e6, P0.L * 1e3
T_AMB = 22.0
F_T, F_N = [0, 1, 2, 3, 4], [0, 400, 600, 857, 0]      # force history inferred from ANSYS average stress


def sigma_of_time(t):                                    # MPa, 0 after step 4
    return float(np.interp(t, F_T, F_N)) / A_MM2


def load_ansys():
    return pd.read_csv(ROOT / "data" / "ansys_3d_wire_merged.csv")


def meff_history(d):
    """1D MEFF driven by the ANSYS T(t) and the inferred force history (steps 1-6)."""
    u = UniaxialMEFF()
    xi, y = [], []
    for t, T in zip(d.time, d["T"]):
        s = sigma_of_time(t)
        u.step(s, T)
        xi.append(u.xi); y.append(L_MM * u.strain(s))
    return np.array(xi), np.array(y)


def plant_prescribed_T(T_path, y_stored):
    """plant recovery map with the ANSYS temperature as input (isolates W1-W3)"""
    m, r, out = RecoveryMap(SMAParams(full_stroke_mm=y_stored)), 0.0, []
    for T in T_path:
        r = max(r, m.r_of_T(T)); out.append(y_stored * (1.0 - r))
    return np.array(out)


def heat_cool(w, y_stored, t_heat=40.0, t_cool=900.0, Tn=0.1, seed=1, ctl_cfg=None):
    """controller heating to full stroke, then I = 0 cooling. Returns logs."""
    cfg = ctl_cfg or CtrlCfg()
    ctl = Controller(w.p, cfg)
    rng = np.random.default_rng(seed + 100)
    E0 = w.E_joule
    h = []
    for k in range(int(t_heat / cfg.Ts)):
        I = ctl.step(1.0, w.T + Tn * rng.standard_normal(), w.measure())
        w.step(I, cfg.Ts, T_inf=T_AMB)
        h.append((k * cfg.Ts, w.T, y_stored - w.stroke_mm(), I))
    E_heat = w.E_joule - E0
    c, dt = [], 0.05
    for k in range(int(t_cool / dt)):
        w.step(0.0, dt, T_inf=T_AMB)
        c.append((t_heat + k * dt, w.T, y_stored - w.stroke_mm()))
    return np.array(h), np.array(c), E_heat


def T_at_fraction(T, dd, y0, frac):
    """first temperature where displacement falls to frac*y0 (linear interpolation)"""
    i = int(np.argmax(dd <= frac * y0))
    if dd[i] > frac * y0:
        return np.nan
    a, b = i - 1, i
    return T[a] + (T[b] - T[a]) * (frac * y0 - dd[a]) / (dd[b] - dd[a])


def run():
    d = load_ansys()
    xi_1d, y_1d = meff_history(d)
    d["xi_1d"], d["y_1d"] = xi_1d, y_1d

    s4 = d[d.step == 4].iloc[-1]
    i4 = s4.name
    xi_st = xi_1d[i4]
    y_st = L_MM * UniaxialMEFF().eps_tr * xi_st            # stored displacement at zero stress

    d5 = d[d.step == 5].reset_index(drop=True)
    y_pl_pres = plant_prescribed_T(d5["T"].values, y_st)

    # ---- cycle 1: controller heating + cooling -----------------------------------------
    p1 = SMAParams(full_stroke_mm=y_st)
    w = SMAWire(p1, T0=35.0, seed=1)
    h1, c1, E1 = heat_cool(w, y_st)
    cool1_end = w.T

    # ---- cycle 2: reload with MEFF at 22 C, then the same heating from 22 C ----------------
    u = UniaxialMEFF(); u.xi = 0.0
    for t in np.linspace(0, 4, 4001):
        s = sigma_of_time(t); u.step(s, T_AMB)
    xi2 = u.xi
    y_st2 = L_MM * u.eps_tr * xi2
    w.r = 1.0 - xi2 / xi_st                                 # hand-off: recovered fraction after reload
    p2 = SMAParams(full_stroke_mm=y_st2)
    w2 = SMAWire(p2, T0=T_AMB, seed=2)
    w2.r = 0.0
    h2, c2, E2 = heat_cool(w2, y_st2, seed=2)

    # ---- metrics ----------------------------------------------------------------------------
    rows = []
    th1 = h1[:, 1]
    for frac, lab in ((0.9, "10 % recovered"), (0.5, "50 % recovered"), (0.1, "90 % recovered")):
        rows.append((lab,
                     T_at_fraction(d5["T"].values, d5.defl_max.values, d5.defl_max.values[0], frac),
                     T_at_fraction(d5["T"].values, d5.y_1d.values, d5.y_1d.values[0], frac),
                     T_at_fraction(d5["T"].values, y_pl_pres, y_st, frac),
                     T_at_fraction(th1, h1[:, 2], y_st, frac)))
    return dict(d=d, d5=d5, y_st=y_st, xi_st=xi_st, y_pl_pres=y_pl_pres, h1=h1, c1=c1, E1=E1, cool1_end=cool1_end,
                y_st2=y_st2, xi2=xi2, h2=h2, c2=c2, E2=E2, rows=rows, s4=s4)


def report(R):
    d, d5, s4 = R["d"], R["d5"], R["s4"]
    A = lambda t: d.iloc[(d.time - t).abs().argmin()]
    print("== 1. loading and unloading (steps 1-4), 1D MEFF vs ANSYS ==")
    for t, lab in ((3.0, "peak 857 N"), (4.0, "end of unload/warm (35 C)")):
        r = A(t)
        print(f"{lab:28s} ANSYS defl_max {r.defl_max:7.3f} mm | 1D {r.y_1d:7.3f} mm | {100 * (r.y_1d / r.defl_max - 1):+5.1f} %")
    print(f"ANSYS average stress still {s4.s_avg:.1f} MPa at end of step 4 (1D: 0); elastic part {s4.s_avg / 46e3 * L_MM:.3f} mm")
    print(f"stored displacement handed to the plant: {R['y_st']:.3f} mm (xi = {R['xi_st']:.3f}); plant default y_max 6.64 mm; ANSYS 6.629 mm")

    print("\n== 2. recovery, step 5: temperature at which the displacement has fallen to the stated level ==")
    print(f"{'':16s} {'ANSYS':>8s} {'1D MEFF':>9s} {'plant,T=ANSYS':>14s} {'plant+ctrl':>11s}   (deg C)")
    for lab, a, b, c, e in R["rows"]:
        print(f"{lab:16s} {a:8.2f} {b:9.2f} {c:14.2f} {e:11.2f}")
    print(f"final displacement after step 5: ANSYS {d5.defl_max.values[-1]:.4f} mm, 1D {d5.y_1d.values[-1]:.4f} mm, "
          f"plant(T=ANSYS) {R['y_pl_pres'][-1]:.4f} mm, plant+controller {R['h1'][-1, 2]:.4f} mm")

    d6 = d[d.step == 6]
    print("\n== 3. cooling, step 6 ==")
    print(f"ANSYS displacement change during cooling: {d6.defl_max.max() - d6.defl_max.min():.2e} mm")
    c1 = R["c1"]
    print(f"plant displacement change during cooling {R['h1'][-1, 2]:.4f} -> {c1[-1, 2]:.4f} mm, "
          f"T {R['h1'][-1, 1]:.1f} -> {R['cool1_end']:.2f} C in {c1[-1, 0] - c1[0, 0]:.0f} s")

    h1 = R["h1"]
    i = np.argmax(h1[:, 2] <= 0.02 * R["y_st"])
    print("\n== 4. controller heating (cycle 1, start 35 C, ambient 22 C) ==")
    print(f"time to 98 % stroke {h1[i, 0]:.2f} s, Tmax {h1[:, 1].max():.2f} C, Imax {h1[:, 3].max():.1f} A, "
          f"electrical energy {R['E1']:.0f} J")

    h2 = R["h2"]; j = np.argmax(h2[:, 2] <= 0.02 * R["y_st2"])
    print("\n== 5. second cycle (reload by MEFF at 22 C, controller heating from 22 C) ==")
    print(f"stored displacement cycle 2 {R['y_st2']:.4f} mm vs cycle 1 {R['y_st']:.4f} mm (identical by construction: MEFF has no cycling drift)")
    print(f"time to 98 % stroke {h2[j, 0]:.2f} s, Tmax {h2[:, 1].max():.2f} C, energy {R['E2']:.0f} J (cycle 1: {R['E1']:.0f} J, from 35 C)")
    print(f"final displacement cycle 2 {h2[-1, 2]:.4f} mm")


def figure(R, outdir):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BG, INK, MUT, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
    B, O, G = "#2a78d6", "#eb6834", "#1baf7a"
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": GRID, "axes.labelcolor": MUT,
                         "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5})
    d, d5, h1, c1 = R["d"], R["d5"], R["h1"], R["c1"]
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 8))
    a = ax[0, 0]
    m = d[d.step <= 4]
    a.plot(m.time, m.defl_max, color=O, lw=2, label="ANSYS")
    a.plot(m.time, m.y_1d, color=B, lw=2, label="1D MEFF")
    a.set_xlabel("ANSYS time [s]"); a.set_ylabel("displacement [mm]"); a.legend(frameon=False)
    a.set_title("1  Load / unload + warm (time axis uses inferred force, timing not comparable)", loc="left", fontsize=10)
    a = ax[0, 1]
    a.plot(d5["T"], d5.defl_max, color=O, lw=2.4, label="ANSYS")
    a.plot(d5["T"], d5.y_1d, color=B, lw=1.6, ls="--", label="1D MEFF (T from ANSYS)")
    a.plot(d5["T"], R["y_pl_pres"], color=G, lw=1.4, ls=":", label="plant map (T from ANSYS)")
    a.plot(h1[:, 1], h1[:, 2], color=INK, lw=1.2, label="plant + controller")
    a.set_xlim(46, 54); a.set_xlabel("temperature [°C]"); a.set_ylabel("displacement [mm]"); a.legend(frameon=False, fontsize=8.5)
    a.set_title("2  Recovery, 35 → 52 °C (zoom 46–54 °C)", loc="left", fontsize=10)
    a = ax[1, 0]
    a.plot(h1[:, 0], h1[:, 2], color=INK, lw=2, label="displacement")
    a.set_xlabel("time [s]"); a.set_ylabel("displacement [mm]")
    b = a.twinx(); b.plot(h1[:, 0], h1[:, 1], color=O, lw=1.6, label="T"); b.set_ylabel("T [°C]", color=O); b.grid(False)
    a.set_title("3  Controller heating (cycle 1)", loc="left", fontsize=10)
    a = ax[1, 1]
    a.plot(c1[:, 0], c1[:, 2], color=INK, lw=2, label="displacement")
    a.set_xlabel("time [s]"); a.set_ylabel("displacement [mm]"); a.set_ylim(-0.1, 0.2)
    b = a.twinx(); b.plot(c1[:, 0], c1[:, 1], color=O, lw=1.6); b.set_ylabel("T [°C]", color=O); b.grid(False)
    a.set_title("4  Cooling, I = 0: T falls, displacement constant (0.047 mm left, see text)", loc="left", fontsize=10)
    fig.tight_layout(); Path(outdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(outdir) / "fig7_full_cycle.png", dpi=160, facecolor=BG)


if __name__ == "__main__":
    R = run()
    report(R)
    if len(sys.argv) > 1:
        figure(R, sys.argv[1])
