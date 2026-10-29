"""1D MEFF loading model against the ANSYS standalone-wire run (steps 1-4, T = 22 C, 857 N).
python scripts/loading_study.py [figure_dir]"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sma_control import SMAParams
from sma_control.loading import UniaxialMEFF

ROOT = Path(__file__).resolve().parents[1]
p = SMAParams()
A, L = p.A * 1e6, p.L * 1e3                      # mm2, mm
XI_SAT = np.sqrt(2 / 3) * p.meff.eps_L           # 0.03266, ANSYS equivalent plastic strain at saturation

d = pd.read_csv(ROOT / "data" / "ansys_3d_wire_merged.csv")
d = d[d.step <= 4].reset_index(drop=True)
d["xi_avg"] = d.p_avg / XI_SAT
d["xi_max"] = d.p_max / XI_SAT

# force history inferred from the average stress (about 400 N at 1 s, 600 N at 2 s, 857 N at 3 s, unloaded at 4 s)
F = np.interp(d.time, [0, 1, 2, 3, 4], [0, 400, 600, 857, 0])
sig = F / A

u = UniaxialMEFF()
xi, eps = [], []
for s, T in zip(sig, d["T"]):
    u.step(s, T)
    xi.append(u.xi); eps.append(u.strain(s))
d["F"], d["sig"], d["xi_1d"], d["y_1d"] = F, sig, xi, np.array(eps) * L


def cross(col, level):                           # ANSYS average stress where xi_avg first reaches `level`
    m = d[(d.step < 4) & (d.time > 1.9)]
    i = np.argmax(m[col].values >= level)
    a, b = m.iloc[i - 1], m.iloc[i]
    return a.s_avg + (b.s_avg - a.s_avg) * (level - a[col]) / (b[col] - a[col])

if __name__ == "__main__":
    print("== transformation stresses, T = 22 C ==")
    print(f"                       1D theory   ANSYS (average stress s_avg)")
    print(f"xi = 0.05  (onset)   {u.sigma_start(22) + 0.05 * 7.1:9.1f}   {cross('xi_avg', 0.05):9.1f}")
    print(f"xi = 0.50            {u.sigma_start(22) + 0.50 * 7.1:9.1f}   {cross('xi_avg', 0.50):9.1f}")
    print(f"xi = 0.80            {u.sigma_start(22) + 0.80 * 7.1:9.1f}   {cross('xi_avg', 0.80):9.1f}")
    print(f"xi = 1 (finish)      {u.sigma_finish(22):9.1f}   (xi_avg reaches 0.94 at 275 MPa and does not reach 1)")
    mx = d[d.step < 4]
    i = np.argmax(mx.p_max.values > 1e-4)
    print(f"local onset (max node): s_max = {mx.s_max.values[i]:.1f} MPa at t = {mx.time.values[i]:.2f} s, "
          f"s_avg = {mx.s_avg.values[i]:.1f} MPa -> clamp stress concentration {mx.s_max.values[i] / mx.s_avg.values[i]:.2f}")
    print("\n== displacement of the free end: ANSYS defl_max vs 1D L*eps ==")
    print(" t [s]   F [N]   ANSYS [mm]   1D [mm]   diff")
    for t in (0.565, 1.0917, 2.0, 2.4259, 2.9059, 3.0, 3.3650, 3.6150, 4.0):
        r = d.iloc[(d.time - t).abs().argmin()]
        print(f"{r.time:6.3f} {r.F:7.0f} {r.defl_max:11.3f} {r.y_1d:9.3f} {100 * (r.y_1d / r.defl_max - 1):+6.1f} %")
    last = d.iloc[-1]
    print(f"\nafter unloading: ANSYS defl_max = {last.defl_max:.3f} mm, 1D stored stroke = {last.y_1d:.3f} mm "
          f"(theory xi = 1: {L * XI_SAT:.3f} mm)")
    print(f"average transformed fraction at end of loading (EPPL avg / saturation) = {d.xi_avg[d.time <= 3].iloc[-1]:.3f}")
    print("\n== stress shift of the recovery window ==")
    for s in (0, 50, 100, 200, 272.8):
        print(f"sigma = {s:6.1f} MPa -> A_s = {u.A_s(s):6.1f} C, A_f = {u.A_f(s):6.1f} C   (window {u.A_f(s) - u.A_s(s):.3f} K)")


def figure(outdir):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BG, INK, MUT, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
    B, O, G = "#2a78d6", "#eb6834", "#1baf7a"
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": GRID, "axes.labelcolor": MUT,
                         "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    m = d[(d.step <= 3) & (d.time > 1.8)]
    s = np.linspace(195, 220, 300)
    th = np.clip([(np.sqrt(2 / 3) * x - p.meff.R / (1 + u.kappa * p.meff.m)) / (p.meff.h * p.meff.eps_L) for x in s], 0, 1)
    ax[0].plot(s, th, color=B, lw=2.5, label="1D MEFF theory")
    ax[0].plot(m.s_avg, m.xi_avg, color=O, lw=2, label="ANSYS (average over wire)")
    ax[0].set_xlim(195, 225); ax[0].set_ylim(-0.02, 1.05)
    ax[0].set_xlabel("average stress [MPa]"); ax[0].set_ylabel("transformed fraction ξ")
    ax[0].legend(frameon=False, loc="center left", fontsize=8.5)
    ax[0].set_title("Transformation at 22 °C: onset 204 → finish 211 MPa", loc="left", fontsize=10)
    ax[1].plot(d.time, d.defl_max, color=O, lw=2, label="ANSYS (free-end node)")
    ax[1].plot(d.time, d.y_1d, color=B, lw=2.2, label="1D MEFF, L·ε")
    ax[1].plot(d.time, d.y_1d * 1.03, color=B, lw=1.2, ls=":", label="1D × 1.03")
    ax[1].set_xlabel("time [s]"); ax[1].set_ylabel("end displacement [mm]"); ax[1].legend(frameon=False, loc="upper left", fontsize=8.5)
    ax[1].set_title("Load 0 → 857 N (0–3 s), unload (3–4 s)", loc="left", fontsize=10)
    sg = np.linspace(0, 280, 100)
    ax[2].fill_between(sg, [u.A_s(x) for x in sg], [u.A_f(x) for x in sg], color=B, alpha=.35)
    ax[2].plot(sg, [u.A_s(x) for x in sg], color=B, lw=1.8)
    ax[2].set_xlabel("constant stress during heating [MPa]"); ax[2].set_ylabel("recovery temperature [°C]")
    ax[2].text(8, 82, "slope 0.146 K/MPa\nwindow stays 1.04 K", color=INK, fontsize=9)
    ax[2].set_title("Stress moves the recovery window", loc="left", fontsize=10)
    fig.tight_layout(); Path(outdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(outdir) / "fig6_loading_vs_ansys.png", dpi=160, facecolor=BG)


if __name__ == "__main__" and len(sys.argv) > 1:
    figure(sys.argv[1])
