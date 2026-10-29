"""Four x-y charts of one full cycle: stress-strain, strain-temperature, current-temperature, current-strain.
Phases: 1 load 22 C (isothermal) | 2 unload | 3 heat with the position controller (controller.py, r_ref = 1) | 4 cool, I = 0.
python scripts/cycle_xy_charts.py [figure_dir]"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sma_control import SMAParams
from sma_control.cycle_plant import SMAWireCycle
from sma_control.controller import Controller, CtrlCfg

T_AMB, F_MAX = 22.0, 857.0
P0 = SMAParams()
SIG_MAX = F_MAX / (P0.A * 1e6)
Y_ST = P0.L * 1e3 * np.sqrt(2 / 3) * P0.meff.eps_L           # stored displacement at xi = 1


def run(t_heat=12.0, t_cool=600.0, seed=1):
    p = SMAParams(full_stroke_mm=Y_ST)
    w = SMAWireCycle(p, T0=T_AMB, xi0=0.0)
    rows = []                                                 # phase, sigma, strain %, T, I, t
    def add(ph, sig, I, t):
        rows.append((ph, sig, 100 * w.u.strain(sig, w.xi), w.T, I, t))
    n = 400
    for k in range(n + 1):                                    # 1 load
        s = SIG_MAX * k / n; w.step_isothermal(s, T_AMB); add(1, s, 0.0, 0.0)
    for k in range(n + 1):                                    # 2 unload
        s = SIG_MAX * (1 - k / n); w.step_isothermal(s, T_AMB); add(2, s, 0.0, 0.0)
    ctl, cfg = Controller(p, CtrlCfg()), CtrlCfg()
    rng = np.random.default_rng(seed)
    for k in range(int(t_heat / cfg.Ts)):                     # 3 heat under the controller
        contraction = Y_ST - w.disp_mm(0.0)
        y_meas = contraction + p.noise_mm * rng.standard_normal()
        I = ctl.step(1.0, w.T + 0.1 * rng.standard_normal(), y_meas)
        w.step(I, cfg.Ts, 0.0, T_AMB); add(3, 0.0, I, (k + 1) * cfg.Ts)
    E_heat = w.E_joule
    dt = 0.05
    for k in range(int(t_cool / dt)):                         # 4 cool
        w.step(0.0, dt, 0.0, T_AMB)
        if k % 10 == 0:
            add(4, 0.0, 0.0, t_heat + (k + 1) * dt)
    return np.array(rows), E_heat


def report(R, E):
    h = R[R[:, 0] == 3]
    i99 = int(np.argmax(h[:, 2] <= 0.01 * 100 * 0.03266))
    print(f"controller heating: Imax {h[:, 4].max():.1f} A, I at 3 s / 4.5 s / 8 s = "
          f"{np.interp(3, h[:, 5], h[:, 4]):.1f} / {np.interp(4.5, h[:, 5], h[:, 4]):.1f} / {np.interp(8, h[:, 5], h[:, 4]):.2f} A")
    print(f"Tmax {h[:, 3].max():.2f} C, final strain {h[-1, 2]:.3f} % (stored 3.266 %), electrical energy over 12 s: {E:.0f} J")
    print(f"strain at peak load {R[R[:, 0] == 1][-1, 2]:.3f} %, after unload {R[R[:, 0] == 2][-1, 2]:.3f} %")


def figure(R, outdir):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BG, INK, MUT, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
    col = {1: "#2a78d6", 2: "#52514e", 3: "#eb6834", 4: "#1baf7a"}
    lab = {1: "1 load (22 °C)", 2: "2 unload (22 °C)", 3: "3 heat, controller", 4: "4 cool, I = 0"}
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": GRID, "axes.labelcolor": MUT,
                         "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 10})
    ph = R[:, 0]
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    def draw(a, xi, yi, phases, xl, yl, title):
        for k in phases:
            m = ph == k
            a.plot(R[m, xi], R[m, yi], color=col[k], lw=2.4, label=lab[k])
        a.set_xlabel(xl); a.set_ylabel(yl); a.set_title(title, loc="left", fontsize=10.5)
    draw(ax[0, 0], 2, 1, (1, 2, 3, 4), "strain ε [%]", "stress σ [MPa]", "Stress against strain")
    ax[0, 0].legend(frameon=False, fontsize=8.5, loc="center", bbox_to_anchor=(0.45, 0.32))
    draw(ax[0, 1], 3, 2, (1, 2, 3, 4), "temperature T [°C]", "strain ε [%]", "Strain against temperature")
    ax[0, 1].annotate("load / unload at 22 °C", (22.3, 1.6), (26, 1.9), color=MUT, fontsize=8.5, arrowprops=dict(arrowstyle="-", color=MUT))
    draw(ax[1, 0], 3, 4, (3, 4), "temperature T [°C]", "current I [A]", "Current against temperature (I = 0 in phases 1, 2, 4)")
    draw(ax[1, 1], 2, 4, (3, 4), "strain ε [%]", "current I [A]", "Current against strain (I = 0 in phases 1, 2, 4)")
    ax[1, 0].legend(frameon=False, fontsize=8.5, loc="center left", bbox_to_anchor=(0.0, 0.45)); ax[1, 1].legend(frameon=False, fontsize=8.5, loc="upper left")
    ax[1, 0].annotate("preheat at the 35 A limit", (38, 35), (30, 31.5), color=MUT, fontsize=8.5)
    ax[1, 0].annotate("recovery window:\ncurrent falls to about 22 A", (50.9, 25), (38, 14), color=MUT, fontsize=8.5, arrowprops=dict(arrowstyle="-", color=MUT))
    ax[1, 0].annotate("hold, then switch off", (51.9, 4.5), (41, 3.5), color=MUT, fontsize=8.5, arrowprops=dict(arrowstyle="-", color=MUT))
    ax[1, 1].annotate("preheat: strain still 3.27 %,\ncurrent rises to 35 A", (3.25, 30), (1.7, 31), color=MUT, fontsize=8.5, arrowprops=dict(arrowstyle="-", color=MUT))
    ax[1, 1].annotate("recovery: strain falls\nat about 22 to 26 A (time runs right to left)", (1.5, 23.5), (0.5, 9), color=MUT, fontsize=8.5, arrowprops=dict(arrowstyle="-", color=MUT))
    fig.tight_layout(); Path(outdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(outdir) / "fig11_cycle_xy.png", dpi=150, facecolor=BG)


if __name__ == "__main__":
    R, E = run(); report(R, E)
    if len(sys.argv) > 1:
        figure(R, sys.argv[1])
