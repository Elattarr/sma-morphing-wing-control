"""Stress, strain, temperature, current and martensite fraction through one full one-way cycle (and the reload).
python scripts/cycle_charts.py [figure_dir]

P1 load 0 -> 857 N at 22 C (isothermal, quasi-static)      P4 cool, I = 0, sigma = 0
P2 unload at 22 C                                           P5 reload at 22 C (isothermal)
P3 heat with 15 A to A_f + 3 K, sigma = 0                   P6 unload at 22 C
The x axis gives every phase the same width (the durations differ by a factor 100); real durations are in the tick labels.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sma_control import SMAParams
from sma_control.cycle_plant import SMAWireCycle

T_AMB, I_HEAT, F_MAX = 22.0, 15.0, 857.0
P = SMAParams()
SIG_MAX = F_MAX / (P.A * 1e6)                                   # 272.8 MPa


def run():
    w = SMAWireCycle(P, T0=T_AMB, xi0=0.0)
    rows, edges = [], []                                          # rows: phase, frac, t, sigma, strain %, T, I, xi
    def log(ph, frac, t, sig, I):
        rows.append((ph, frac, t, sig, 100 * w.u.strain(sig, w.xi), w.T, I, w.xi))

    def isothermal(ph, s0, s1, dur, n=600):
        for k in range(n + 1):
            sig = s0 + (s1 - s0) * k / n
            w.step_isothermal(sig, T_AMB)
            log(ph, k / n, dur * k / n, sig, 0.0)

    isothermal(1, 0.0, SIG_MAX, 60.0); isothermal(2, SIG_MAX, 0.0, 60.0)
    # P3 heating
    Toff, t, dt, k = w.u.A_f(0.0) + 3.0, 0.0, 0.01, 0
    t_list = []
    while w.T < Toff:
        w.step(I_HEAT, dt, 0.0, T_AMB); t += dt; t_list.append((t, w.T, w.xi))
    t_heat = t
    for t_, T_, xi_ in t_list:                                    # replay for logging at uniform fractions
        pass
    w2 = SMAWireCycle(P, T0=T_AMB, xi0=1.0)
    t = 0.0
    while w2.T < Toff:
        w2.step(I_HEAT, dt, 0.0, T_AMB); t += dt
        rows.append((3, t / t_heat, t, 0.0, 100 * w2.u.strain(0.0, w2.xi), w2.T, I_HEAT, w2.xi))
    w = w2
    t_cool, dtc = 600.0, 0.05
    for k in range(int(t_cool / dtc)):
        w.step(0.0, dtc, 0.0, T_AMB)
        if k % 10 == 0:
            rows.append((4, (k + 1) * dtc / t_cool, (k + 1) * dtc, 0.0, 100 * w.u.strain(0.0, w.xi), w.T, 0.0, w.xi))
    w.T = T_AMB                                                   # next phases are isothermal at ambient (T has relaxed to 22.0x)
    isothermal(5, 0.0, SIG_MAX, 60.0); isothermal(6, SIG_MAX, 0.0, 60.0)
    return np.array(rows), t_heat


def report(R, t_heat):
    def end(ph):
        r = R[R[:, 0] == ph][-1]
        return r
    names = {1: "after load", 2: "after unload", 3: "after heating", 4: "after cooling", 5: "after reload", 6: "after 2nd unload"}
    print(f"{'phase end':18s} {'sigma [MPa]':>11s} {'strain [%]':>10s} {'T [C]':>7s} {'xi':>6s}")
    for ph in range(1, 7):
        r = end(ph)
        print(f"{names[ph]:18s} {r[3]:11.1f} {r[4]:10.3f} {r[5]:7.2f} {r[7]:6.3f}")
    h = R[R[:, 0] == 3]
    print(f"\nheating at {I_HEAT:.0f} A takes {t_heat:.1f} s, Tmax {h[:, 5].max():.1f} C; cooling phase 600 s")
    print(f"peak strain at {SIG_MAX:.1f} MPa: {R[R[:, 0] == 1][-1, 4]:.3f} % ; stored strain after unload: {R[R[:, 0] == 2][-1, 4]:.3f} %")


def figure(R, t_heat, outdir):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BG, INK, MUT, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
    C = dict(sig="#2a78d6", eps="#1baf7a", T="#eb6834", I="#7a4fd1", xi=INK)
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": GRID, "axes.labelcolor": MUT,
                         "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5})
    x = R[:, 0] - 1 + R[:, 1]
    names = ["1 load\n60 s, 22 °C", "2 unload\n60 s, 22 °C", f"3 heat, 15 A\n{t_heat:.1f} s", "4 cool, I = 0\n600 s",
             "5 reload\n60 s, 22 °C", "6 unload\n60 s, 22 °C"]
    fig = plt.figure(figsize=(15, 10.5))
    gl = fig.add_gridspec(5, 1, left=0.06, right=0.60, top=0.95, bottom=0.08, hspace=0.12)
    gr = fig.add_gridspec(2, 1, left=0.69, right=0.97, top=0.95, bottom=0.08, hspace=0.32)
    series = [("stress σ [MPa]", 3, C["sig"]), ("strain ε [%]", 4, C["eps"]), ("temperature T [°C]", 5, C["T"]),
              ("current I [A]", 6, C["I"]), ("martensite fraction ξ", 7, C["xi"])]
    axs = []
    for i, (lab, col, c) in enumerate(series):
        a = fig.add_subplot(gl[i, 0], sharex=axs[0] if axs else None)
        a.plot(x, R[:, col], color=c, lw=2)
        a.set_ylabel(lab)
        for k in range(1, 6):
            a.axvline(k, color=GRID, lw=1.2)
        axs.append(a)
        if i < 4:
            a.tick_params(labelbottom=False)
    axs[-1].set_xticks(np.arange(6) + 0.5); axs[-1].set_xticklabels(names, fontsize=8.5)
    axs[0].set_title("One full cycle: load, unload, heat, cool, reload (phases drawn equally wide, real durations below)", loc="left", fontsize=10)
    ph = R[:, 0]
    a = fig.add_subplot(gr[0])
    for k, col in ((1, C["sig"]), (2, "#52514e"), (3, C["T"])):
        m = ph == k; a.plot(R[m, 4], R[m, 3], color=col, lw=2.2)
    a.annotate("1 load: elastic, then plateau at 204–211 MPa", (0.45, 205), (0.2, 235), color=C["sig"], fontsize=8.5)
    a.annotate("2 unload (elastic)", (3.5, 120), (2.2, 90), color="#52514e", fontsize=8.5)
    a.annotate("3 heat at σ = 0: strain falls to 0", (1.6, 2), (0.3, 30), color=C["T"], fontsize=8.5)
    a.set_xlabel("strain ε [%]"); a.set_ylabel("stress σ [MPa]"); a.set_title("Stress against strain", loc="left", fontsize=10)
    a = fig.add_subplot(gr[1])
    for k, lab in ((3, "3 heating"), (4, "4 cooling")):
        m = ph == k; a.plot(R[m, 5], R[m, 4], color={3: C["T"], 4: "#2a78d6"}[k], lw=2.2, label=lab)
    a.set_xlabel("temperature T [°C]"); a.set_ylabel("strain ε [%]"); a.legend(frameon=False, loc="center left")
    a.set_title("Strain against temperature at σ = 0", loc="left", fontsize=10)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(outdir) / "fig10_cycle_charts.png", dpi=150, facecolor=BG, bbox_inches="tight")


if __name__ == "__main__":
    R, th = run(); report(R, th)
    if len(sys.argv) > 1:
        figure(R, th, sys.argv[1])
