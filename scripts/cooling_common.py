"""Shared helpers for cooling_one_way.py and cooling_two_way.py"""
import numpy as np

T_AMB, I_HEAT = 22.0, 15.0
BG, INK, MUT, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"


def cycle(w, sigma, t_cool=900.0):
    """heat with 15 A until A_f(sigma)+3 K, then I = 0 for t_cool s. Returns log (t, T, xi, d, heating?) and index of end of heating."""
    Toff = w.u.A_f(sigma) + 3.0
    log, t, dt = [], 0.0, 0.02
    while True:
        w.step(I_HEAT if w.T < Toff else 0.0, dt, sigma, T_AMB); t += dt
        log.append((t, w.T, w.xi, w.disp_mm(sigma), 1))
        if w.T >= Toff or t > 600:
            break
    nh, dt = len(log), 0.05
    for _ in range(int(t_cool / dt)):
        w.step(0.0, dt, sigma, T_AMB); t += dt
        log.append((t, w.T, w.xi, w.disp_mm(sigma), 0))
    return np.array(log), nh


def first_forward(log, nh):
    """temperature where xi first rises 0.001 above its minimum during cooling; nan if it never does"""
    c = log[nh:]
    xmin = c[:, 2].min()
    return c[int(np.argmax(c[:, 2] > xmin + 1e-3)), 1] if c[-1, 2] > xmin + 1e-3 else np.nan


def style():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": GRID, "axes.labelcolor": MUT,
                         "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5})
    return plt


def time_panel(ax, L, title):
    ax.plot(L[:, 0], L[:, 3], color=INK, lw=2)
    ax.set_xlabel("time [s]"); ax.set_ylabel("displacement [mm]")
    b = ax.twinx(); b.plot(L[:, 0], L[:, 1], color="#eb6834", lw=1.6)
    b.set_ylabel("T [°C]", color="#eb6834"); b.grid(False)
    ax.set_title(title, loc="left", fontsize=10)
