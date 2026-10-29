"""Parameter identification on simulated bench data: true wire (parameters differ from the nominal ones) -> tests A, B, C -> fitted values.

python scripts/identify_bench.py        (0.5 mm wire, 20 N/mm... see below: spring 1.25 N/mm, 250 MPa preload, 3 A limit)
Sensors: thermocouple noise 0.1 K and first-order lag, current noise 0.3 %, voltage noise 0.5 %, displacement noise 0.02 mm, 20 Hz.
"""
import os
import sys
from dataclasses import replace
from pathlib import Path
import numpy as np
os.environ.setdefault("SMA_D_MM", "0.5"); os.environ.setdefault("SMA_IMAX", "3"); os.environ.setdefault("SMA_SPRING_K", "1.25")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_controller as ev
from sma_control.cycle_controller import HomingRoutine
from sma_control.identify import fit_thermal, fit_latent, fit_resistance_contracted

NOM = ev.NOM
TS = 0.05


class Bench:
    """measurement chain: lagged thermocouple, noisy current and voltage, noisy displacement"""
    def __init__(self, w, tau_lag, seed):
        self.w, self.tau, self.rng, self.Tf = w, tau_lag, np.random.default_rng(seed), w.T

    def measure(self, I_applied):
        w, rng = self.w, self.rng
        self.Tf += TS / (self.tau + TS) * (w.T - self.Tf)
        Tm = self.Tf + 0.1 * rng.standard_normal()
        Im = I_applied * (1 + 0.003 * rng.standard_normal())
        V = I_applied * w.p.resistance(1.0 - w.xi) * (1 + 0.005 * rng.standard_normal())
        d = w.disp_mm(ev.SIG) + 0.02 * rng.standard_normal()
        return Tm, Im, V, d


def test_A(p_true, tau_lag, seed):
    w = ev.make_plant(p_true, 22.0, 1.0, ev.SIG)
    b = Bench(w, tau_lag, seed)
    I_test = float(np.sqrt(1.0 * NOM.heat_capacity / NOM.resistance(0.0)))        # about 1 K/s with the NOMINAL numbers
    rec, heating, t_end = [], True, None
    I_app, T_REST = 0.0, 10.0                                   # 10 s at rest first: ambient is measured on the cold, unpowered wire
    for k in range(int(1500 / TS)):
        t = k * TS
        Tm, Im, V, d = b.measure(I_app)
        rec.append((t, Im, V, Tm))
        if t >= T_REST and heating and (Tm - rec[0][3] > 12.0 or t > 300.0 + T_REST):
            heating, t_end = False, t - T_REST
        I_app = I_test if (heating and t >= T_REST) else 0.0
        ev.plant_step(w, I_app, TS, ev.SIG, 22.0)
        if not heating and (Tm - 22.0 < 1.5 or t > 3.0 * t_end + 60.0 + T_REST):
            break
    return np.array(rec)


def test_B_C(p_true, tau_lag, seed):
    """the homing heating record"""
    cfg = ev.make_cfg()
    w = ev.make_plant(p_true, 22.0, 1.0, ev.SIG)
    ctl = ev.make_ctl(ev.SIG, cfg)
    h = HomingRoutine(ctl)
    b = Bench(w, tau_lag, seed)
    I_app, rec = 0.0, []
    while not h.done and h.t < 3000:
        Tm, Im, V, d = b.measure(I_app)
        I_cmd = h.step(Tm, d)
        if h.phase in ("wait", "heat") or h.t < 1.0:
            rec.append((h.t, Im, V, Tm, d))
        I_app = I_cmd
        ev.plant_step(w, I_app, TS, ev.SIG, 22.0)
    return np.array(rec), h.cal


def run_case(name, p_true, tau_lag, seed=3):
    A = test_A(p_true, tau_lag, seed)
    th = fit_thermal(A[:, 0], A[:, 1], A[:, 2], A[:, 3])
    rec, cal = test_B_C(p_true, tau_lag, seed + 1)
    if cal is None:
        return name, None
    L = NOM.L * 1e3
    d_ext, d_con = cal.eps_ext * L, cal.eps_con * L
    # heat phase only: from the first power-on to the end
    on = np.argmax(rec[:, 1] > 0.05)
    t, I, V, Tm, d = (rec[on - 1:, j] - (rec[on - 1, 0] if j == 0 else 0.0) for j in range(5))
    mL = fit_latent(t, I, V, Tm, d, th, d_ext, d_con)
    RA = fit_resistance_contracted(V, I, d, d_ext, d_con)
    return name, dict(mc=(th.mc, p_true.heat_capacity), hS=(th.hS, p_true.h * p_true.surface), tau=(th.tau_lag, tau_lag),
                      RM=(th.R, p_true.resistance(0.0)), RA=(RA, p_true.resistance(1.0)), mL=(mL, p_true.latent_total), rms=th.rms)


CASES = [("nominal, lag 0.5 s", NOM, 0.5), ("c = 320 (m c x0.38)", replace(NOM, c=320.0), 0.5), ("h = 5", replace(NOM, h=5.0), 0.5),
         ("h = 82", replace(NOM, h=82.0), 0.5), ("latent x0.3", replace(NOM, latent_J_per_mm3=0.3 * NOM.latent_J_per_mm3), 0.5),
         ("latent x2", replace(NOM, latent_J_per_mm3=2.0 * NOM.latent_J_per_mm3), 0.5),
         ("resistivity x1.2", replace(NOM, rho_eM=1.2 * NOM.rho_eM, rho_eA=1.2 * NOM.rho_eA), 0.5), ("lag 0.1 s", NOM, 0.1), ("lag 3 s", NOM, 3.0)]

if __name__ == "__main__":
    print("== PARAMETER IDENTIFICATION (simulated bench data; fitted value vs true value, error in %) ==")
    print(f"{'case':24s} " + " ".join(f"{k:>16s}" for k in ("m c [J/K]", "h S [W/K]", "tau lag [s]", "R_M [ohm]", "R_A [ohm]", "m Lambda [J]")) + "   rms [K]")
    worst = {}
    for name, p, lag in CASES:
        n, r = run_case(name, p, lag)
        if r is None:
            print(f"{name:24s} homing failed"); continue
        cells = []
        for k in ("mc", "hS", "tau", "RM", "RA", "mL"):
            f, t = r[k]; e = 100 * (f - t) / t
            worst[k] = max(worst.get(k, 0.0), abs(e)); cells.append(f"{f:8.4g} {e:+6.1f}%")
        print(f"{name:24s} " + " ".join(f"{c:>16s}" for c in cells) + f"   {r['rms']:.3f}")
    print("worst |error| % : " + ", ".join(f"{k} {v:.1f}" for k, v in worst.items()))
