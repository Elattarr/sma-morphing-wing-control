"""Figures for the parameter sweep: python scripts/sweep_figures.py <outdir>"""
import sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from param_sweep import *

out = Path(sys.argv[1]); out.mkdir(exist_ok=True, parents=True)
BG, INK, MUT, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#e34948"]
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": GRID, "axes.labelcolor": MUT,
                     "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 10})

def trace(p, r_ref, T_amb=22.0, T_end=12.0, seed=1):
    cfg = CtrlCfg(); w = SMAWire(p, T0=T_amb, seed=seed); ctl = Controller(NOM, cfg)
    rng = np.random.default_rng(seed + 100); n = int(T_end / cfg.Ts)
    a = np.zeros((n, 4))
    for k in range(n):
        Im = ctl.step(r_ref, w.T + 0.1 * rng.standard_normal(), w.measure())
        w.step(Im, cfg.Ts, T_inf=T_amb)
        a[k] = (k * cfg.Ts, w.r, w.T, Im)
    return a

# Fig 1: nominal response, both steps
fig, ax = plt.subplots(3, 1, figsize=(7.5, 7), sharex=True)
for rr, col in ((0.5, C[0]), (1.0, C[1])):
    a = trace(NOM, rr)
    ax[0].plot(a[:, 0], a[:, 1] * 6.64, color=col, lw=2, label=f"reference {rr:.1f}")
    ax[0].axhline(rr * 6.64, color=col, lw=1, ls=":")
    ax[1].plot(a[:, 0], a[:, 2], color=col, lw=2)
    ax[2].plot(a[:, 0], a[:, 3], color=col, lw=2)
ax[1].axhspan(NOM.meff.As, NOM.meff.Af, color="#cde2fb", alpha=.8, lw=0); ax[1].text(0.1, 56.5, "shaded band: A_s to A_f (1.04 K)", color=MUT, va="center", fontsize=9)
ax[0].set_ylabel("stroke [mm]"); ax[1].set_ylim(20, 58); ax[1].set_ylabel("wire temperature [°C]"); ax[2].set_ylabel("current [A]"); ax[2].set_xlabel("time [s]")
ax[0].legend(frameon=False, loc="lower right"); ax[0].set_title("Nominal plant: closed-loop step to 50 % and 100 % stroke", loc="left")
fig.tight_layout(); fig.savefig(out / "fig1_nominal_steps.png", dpi=160)

# Fig 2: 50 % step under the three worst perturbations
cases = [("nominal", NOM, C[0]), ("latent heat 46 J", replace(NOM, latent_J_per_mm3=LAM_CC), C[1]),
         ("h_c = 5", replace(NOM, h=5.0), C[2]), ("window 12 K", replace(NOM, meff=wide_window(12.0)), C[3])]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for n_, p, col in cases:
    a = trace(p, 0.5, T_end=10); ax.plot(a[:, 0], a[:, 1], color=col, lw=2, label=n_)
ax.axhline(0.5, color=INK, lw=1, ls=":"); ax.set_ylim(0, 1.02)
ax.set_xlabel("time [s]"); ax.set_ylabel("recovered fraction r"); ax.legend(frameon=False, loc="lower right")
ax.set_title("50 % step: where the nominal controller fails (one-way, overshoot is permanent)", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig(out / "fig2_worst_cases.png", dpi=160)

# Fig 3: overshoot bars, 50 % step, all cases
names, os_ = [], []
for name, p, Ta in CASES:
    names.append(name); os_.append(run(p, 0.5, T_amb=Ta)["os"])
fig, ax = plt.subplots(figsize=(7.5, 5.2))
y = np.arange(len(names))[::-1]
cols = [C[3] if v > 3 else C[0] for v in os_]
ax.barh(y, os_, color=cols, height=.62)
ax.set_yticks(y); ax.set_yticklabels(names); ax.axvline(3, color=INK, lw=1, ls=":"); ax.text(3.6, y[-1] - 0.75, "3 % limit", color=MUT, va="center", fontsize=9); ax.set_ylim(-1.3, len(names) - .4)
for yy, v in zip(y, os_): ax.text(max(v, 3) + 1.2, yy, f"{v:.1f}", va="center", color=INK, fontsize=9)
ax.set_xlabel("overshoot on the 50 % step [% of target]"); ax.grid(axis="y", visible=False)
ax.set_title("Overshoot of the 50 % step across the parameter sweep", loc="left")
fig.tight_layout(); fig.savefig(out / "fig3_overshoot_sweep.png", dpi=160)
print("ok")
