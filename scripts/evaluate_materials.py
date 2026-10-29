"""Other SMA materials without retuning: the same controller, gains and homing on wires whose MEFF constants differ from specimen C.

python scripts/evaluate_materials.py
For every alloy the bias preload is set by the return rule  T_fwd(1) = T_ambient + 12 K  (closed form, see spring.py):
    sigma_ext = ((T_target - T0) beta + R/(1+kappa m) + h eps_L) / sqrt(2/3)
Two controllers per alloy, identical gains:
    'nominal'  built from specimen C constants  (the controller does not know the alloy: homing has to find the windows)
    'own'      built from the alloy's constants (windows from the model, then corrected by homing)
Only the MEFF constants change; density, specific heat, latent heat, resistivities, convection and geometry stay those of the 0.5 mm wire
(a different alloy would also change those: not covered). A T0 shift is equivalent to a stress shift, so the two T0 rows keep the
preload of specimen C (re-solving the preload would cancel the shift exactly). Constants other than those named in the source are [ASSUMED].
"""
import math
import os
import sys
from dataclasses import replace
from pathlib import Path
os.environ.setdefault("SMA_D_MM", "0.5"); os.environ.setdefault("SMA_IMAX", "3"); os.environ.setdefault("SMA_SPRING_K", "1.25")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_controller as ev

M0 = ev.NOM.meff
T_AMB, MARGIN = 22.0, 12.0


def preload(M):
    kappa = math.sqrt(2.0 / 27.0)
    return ((T_AMB + MARGIN - M.T0) * M.beta + M.R / (1.0 + kappa * M.m) + M.h * M.eps_L) / math.sqrt(2.0 / 3.0)


MATERIALS = [
    ("specimen C (baseline)", M0),
    ("Auricchio 2009 type: h 100, eps_L 5.1 %, beta 7.7", replace(M0, h=100.0, eps_L=0.051, beta=7.7)),
    ("cold-worked: h 700 (window 5 K)", replace(M0, h=700.0)),
    ("cold-worked: h 1400 (window 10 K)", replace(M0, h=1400.0)),
    ("T0 +13 K (35 C), preload of specimen C", replace(M0, T0=35.0)),
    ("T0 -10 K (12 C), preload of specimen C", replace(M0, T0=12.0)),
    ("beta 4.5", replace(M0, beta=4.5)),
    ("beta 8.0", replace(M0, beta=8.0)),
    ("eps_L 3 %", replace(M0, eps_L=0.03)),
    ("eps_L 6 %", replace(M0, eps_L=0.06)),
    ("R 120 (narrow hysteresis)", replace(M0, R=120.0)),
    ("R 220 (wide hysteresis)", replace(M0, R=220.0)),
]


FIXED_PRELOAD = ('T0 +13 K', 'T0 -10 K')          # temperature shift at the preload of specimen C: windows move, nothing re-solved


def row(name, p_true, sig, p_ctl):
    L, seg, home = ev.sequence(p_true=p_true, sigma_true=sig, sigma_ctl=sig, p_ctl=p_ctl, label=name)
    s = ev.seq_summary(L, seg, home)
    note = ("FAULT " + home["fault"]) if home["fault"] else ""
    cal = home.get("cal")
    return s, note, cal


if __name__ == "__main__":
    print("== OTHER MATERIALS, gains unchanged: sequence r = 1.0 -> 0.5 -> 0.8 -> 0.0, preload from the return rule ==")
    print(f"{'material':52s} {'preload':>8s} {'A_s':>6s} {'A_f':>6s} | {'nominal ctl':>14s} {'hold err %':>10s} | {'own ctl':>14s} {'hold err %':>10s} {'Tmax':>6s}")
    for name, M in MATERIALS:
        p_true = replace(ev.NOM, meff=M)
        sig = preload(M0) if name.startswith(FIXED_PRELOAD) else preload(M)
        sA, nA, cal = row(name, p_true, sig, None)
        sB, nB, _ = row(name, p_true, sig, p_true)
        As, Af = (cal.As, cal.Af) if cal else (float("nan"),) * 2
        print(f"{name:52s} {sig:7.0f}  {As:6.1f} {Af:6.1f} | {sA['reached']:>6d} of 4 {nA:>0s}{sA['max_err']:10.2f} | {sB['reached']:>6d} of 4 {sB['max_err']:10.2f} {sB['Tmax']:6.1f} {nB}")
    print(f"safety monitor: " + ("no trips in any run" if not ev.TRIPS else "TRIPS " + "; ".join(f"{a[:30]}: {b} at {c:.0f} s" for a, b, c in ev.TRIPS)))
