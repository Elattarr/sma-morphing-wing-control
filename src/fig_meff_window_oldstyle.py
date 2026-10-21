"""Replacement for thesis Fig. 4.5 in the OLD layout, with corrected As/Af from the MEFF closed form.
Curve = ANSYS standalone wire, EPPL max (3.27 % frozen), heating after force removal."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
root = Path(__file__).resolve().parents[1]
m = pd.read_csv(root / "data" / "ansys_3d_wire_merged.csv")
s = m[m.step.isin([4, 5])].sort_values("time")
T, eps = s["T"].values, s["p_max"].values * 100
T0, R, H, beta, eL = 22.0, 167.0, 145.0, 5.6, 0.04
As, Af = T0 + (R - H*eL)/beta, T0 + R/beta
dark, blue = "#990000", "#0033aa"
fig, ax = plt.subplots(figsize=(11, 6.3))
ax.set_facecolor("#fafaff")
ax.axvspan(As, Af, color="#ffd9d9", alpha=.7, zorder=0)
for x in (As, Af): ax.axvline(x, color="red", ls="--", lw=.8, zorder=1)
# loading / unloading at T0, then frozen plateau, recovery from ANSYS, flat at zero
Tp = np.concatenate([[T0, T0], T, [Af + 0.2, 56.0]])
ep = np.concatenate([[0, eps[0]], eps, [0, 0]])
ax.plot(Tp, ep, color=dark, lw=3, label="ANSYS MEFF Simulation", zorder=3)
ax.plot([T0], [eps[0]], "o", color=dark, ms=8, zorder=4)
ax.plot([Af], [0], "D", color=dark, ms=9, zorder=4)
ax.plot([T0], [4.0], "o", color=blue, ms=10)
ax.annotate("$\\varepsilon_L$ = 4% (theoretical max)", (T0, 4.0), xytext=(-2.5, 12), textcoords="offset points", color=blue, fontsize=11, ha="left")
ax.annotate("Frozen martensite   ($\\varepsilon$ = 3.27%)", (T0, eps[0]), xytext=(8, 10), textcoords="offset points", color=blue, fontsize=11)
ax.annotate("Loading / unloading\n$F$ = 857 N,  $T_0$ = 22 °C", (T0, 1.9), xytext=(8, 0), textcoords="offset points", color=blue, fontsize=11)
ax.annotate(f"$A_s$ = {As:.1f} °C", (As, 0), xytext=(-82, 8), textcoords="offset points", color="red", fontsize=11)
ax.annotate(f"$A_f$ = {Af:.1f} °C\nFull recovery", (Af, 0), xytext=(14, 18), textcoords="offset points", color="red", fontsize=11)
ax.set_xlim(19, 59); ax.set_ylim(-0.3, 4.45)
ax.set_xlabel("Temperature [°C]", fontsize=13); ax.set_ylabel("SMA Wire Strain [%]", fontsize=13)
ax.set_title(f"Strain-Temperature Response: Complete SME Cycle\nNiTi SMA Wire — Auricchio-Petrini MEFF Model  ($A_s$ = {As:.1f} °C,  $A_f$ = {Af:.1f} °C,  $T_0$ = 22 °C)", fontsize=13, fontweight="bold")
ax.grid(ls="--", alpha=.35); ax.legend(loc="upper right", fontsize=11)
ins = ax.inset_axes([0.30, 0.12, 0.36, 0.40]); ins.set_facecolor("#fafaff")
ins.axvspan(As, Af, color="#ffd9d9", alpha=.7)
ins.plot(T, eps, color=dark, lw=2); ins.set_xlim(49.5, 52.5); ins.set_ylim(-0.1, 3.5)
for x in (As, Af): ins.axvline(x, color="red", ls="--", lw=.8)
ins.set_title("Zoom on recovery window", fontsize=9); ins.set_xlabel("Temperature [°C]", fontsize=8); ins.tick_params(labelsize=8); ins.grid(ls="--", alpha=.3)
plt.tight_layout(); plt.savefig(root / "figs" / "fig4_5_new_oldstyle.png", dpi=170)
print(As, Af)
