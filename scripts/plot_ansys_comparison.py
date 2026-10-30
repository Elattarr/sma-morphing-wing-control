"""Figure: ANSYS Up 40 mm recovery vs the controller model and the fitted plant, the TE map, and the staircase on the fitted plant.

PYTHONPATH=src:scripts python scripts/plot_ansys_comparison.py [out.png]
Needs docs/ansys_up40_recovery.csv (T, wire-end displacement, TE deflection, wire recovery from the ANSYS run).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import evaluate_heatonly as h            # sets the 2 mm plant and patches the harness (import first)
import evaluate_controller as ev
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sma_control.te_map import TEMap
from sma_control.params import SMAParams
from sma_control.spring import BiasSpring
from sma_control.elements import SMAWireElements

ROOT = Path(__file__).resolve().parents[1]
FIT = dict(spread_K=0.4, tail_frac=0.05, tail_K=5.0)          # fitted to the ANSYS recovery curve (RMS 1.8 % of the stroke)


def quasi_static(P, st, Tg):
    w = SMAWireElements(P, BiasSpring(2.2, 6.36), T0=22.0, xi0=1.0, n=41, **st)
    out = []
    for t in Tg:
        xi = w.xi_next(w.xi, 0.0, t)
        w.xe = w._c(xi, t, w.xe)
        w.xi = float(w.xe.mean())
        out.append(1 - w.xi)
    return np.array(out)


def main(out):
    P = SMAParams(d=2e-3, Imax=35.0)
    a = pd.read_csv(ROOT / "docs/ansys_up40_recovery.csv")
    m = a["T_C"] >= 48
    Tg = np.arange(48, 70, 0.05)
    rf, rs = quasi_static(P, FIT, Tg), quasi_static(P, dict(spread_K=0.0), Tg)
    L, seg, home = ev.sequence(targets=h.UP, hold_s=10.0, t_max=150.0, struct=FIT, label="fit")
    te = TEMap()
    t, r, tg = L[:, 0], L[:, 2], L[:, 3]
    f, ax = plt.subplots(2, 2, figsize=(12, 8)); ax = ax.ravel()
    ax[0].plot(a.T_C[m], a.r_wire_end[m], "k", lw=2, label="ANSYS wire end")
    ax[0].plot(Tg, rs, "--", c="C1", label="controller model (sharp)")
    ax[0].plot(Tg, rf, c="C2", label="fitted plant (0.4 K, 5 % tail 5 K)")
    ax[0].set(xlim=(48, 66), xlabel="T [°C]", ylabel="recovery r", title="Recovery vs temperature"); ax[0].legend(); ax[0].grid(alpha=.3)
    rr = np.linspace(0, 1, 50)
    ax[1].plot(rr, [te.r_to_te(x) for x in rr], "k", lw=2, label="ANSYS map")
    ax[1].plot(rr, -6.632 + rr * (-42.239 + 6.632), "--", c="C1", label="linear")
    ax[1].set(xlabel="wire recovery r", ylabel="TE deflection [mm]", title="TE map (concave)"); ax[1].legend(); ax[1].grid(alpha=.3)
    ax[2].plot(t, [te.r_to_te(x) for x in tg], "--", c="C1", label="target")
    ax[2].plot(t, [te.r_to_te(x) for x in r], "k", lw=2, label="TE deflection (model)")
    ax[2].set(xlabel="t [s]", ylabel="TE [mm]", title="Staircase on ANSYS-fitted plant"); ax[2].legend(); ax[2].grid(alpha=.3)
    ax[3].plot(t, L[:, 1], c="C3"); ax[3].set(xlabel="t [s]"); ax[3].set_ylabel("wire T [°C]", color="C3")
    b = ax[3].twinx(); b.plot(t, L[:, 4], c="C0", lw=.8); b.set_ylabel("current [A]", color="C0")
    ax[3].set_title("Temperature (peak %.1f °C) and current" % L[:, 1].max()); ax[3].grid(alpha=.3)
    plt.tight_layout(); plt.savefig(out, dpi=130)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "docs/figs_wing/heatonly_ansys_comparison.png"))
