"""Preview of the replacement for thesis Fig. 4.5: MEFF recovery window of the standalone wire.
Data: ANSYS standalone 3D wire, heating step (T 22.3 -> 52 C). Lines: closed-form As, Af (Auricchio & Petrini 2004).
  As = T0 + (R - H*eps_L)/beta ,  Af = T0 + R/beta
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
root = Path(__file__).resolve().parents[1]
m = pd.read_csv(root / "data" / "ansys_3d_wire_merged.csv")
s = m[m.step.isin([4, 5])]                       # heating after force removal
T0, R, H, beta, eL = 22.0, 167.0, 145.0, 5.6, 0.04
As, Af = T0 + (R - H*eL)/beta, T0 + R/beta
print(f"As = {As:.2f} C, Af = {Af:.2f} C, window = {Af-As:.2f} K")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for a, xl in zip(ax, [(20, 55), (49, 53)]):
    a.plot(s["T"], s.p_max*100, "--", color="tab:red", lw=1.4, label="ANSYS, EPPL max")
    a.plot(s["T"], s.p_avg*100, "-", color="tab:blue", lw=2, label="ANSYS, EPPL average")
    a.axvspan(As, Af, color="gold", alpha=.35, label=f"MEFF window ({As:.1f}-{Af:.1f} °C)")
    a.axvline(As, color="k", ls=":", lw=1); a.axvline(Af, color="k", ls=":", lw=1)
    a.set_xlim(*xl); a.set_xlabel("Wire temperature [°C]"); a.set_ylabel("Transformation strain EPPL [%]"); a.grid(alpha=.3)
ax[0].legend(loc="lower left", fontsize=8); ax[0].set_title("Heating after force removal", fontsize=10)
ax[1].set_title("Zoom on recovery window", fontsize=10)
ax[1].annotate(f"$A_s$ = {As:.1f}", (As, 3.3), xytext=(-48, 0), textcoords="offset points", fontsize=9)
ax[1].annotate(f"$A_f$ = {Af:.1f}", (Af, 3.3), xytext=(4, 0), textcoords="offset points", fontsize=9)
plt.tight_layout(); plt.savefig(root / "figs" / "fig4_5_new_preview.png", dpi=170)
