"""Compare a BENCH log of the 0.5 mm wire with the model and print pass / fail per criterion.  This is the step that turns "verified" into "validated".

python scripts/validate_bench.py log.csv [--ambient 22]     log.csv: header  t,I,V,T1,d   (s, A, V at the wire, C, mm), >= 20 Hz, one heating run from cold then cooling
python scripts/validate_bench.py --selftest                 synthetic log from a perturbed model (checks the script itself, NOT the wire)

The model is driven by the MEASURED current I(t) (open loop), so the controller is not involved: this tests the plant model only.
Free parameters: thermocouple lag tau (fitted, 0..3 s).  Nothing else is fitted; fit h, m c, R afterwards with scripts/identify_bench.py if a criterion fails.
Criteria (set before the data, change them only with a reason):
  1  RMS temperature error over the run                 <= 3 K
  2  temperature at which the stroke reaches 10 %       within 2 K of the model       (A_s,sigma)
  3  temperature at which the stroke reaches 90 %       within 3 K of the model       (A_f,sigma)
  4  final stroke                                       within 10 % of the model
  5  wire resistance V/I, cold and hot                  within 5 % of the model
"""
import os
import sys
from pathlib import Path
os.environ.setdefault("SMA_D_MM", "0.5"); os.environ.setdefault("SMA_IMAX", "3"); os.environ.setdefault("SMA_SPRING_K", "1.25")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataclasses import replace
import numpy as np
import evaluate_controller as ev

CRIT = dict(rms_T=3.0, dT_start=2.0, dT_finish=3.0, stroke_rel=0.10, R_rel=0.05)


def simulate(t, I, T_amb, p=None, tau=0.0):
    p = p or ev.NOM
    w = ev.make_plant(p, T_amb, 1.0, ev.SIG)
    out, Tf, d0 = [], T_amb, None
    for k in range(len(t)):
        dt = (t[k] - t[k - 1]) if k else 0.05
        Tf += dt / (tau + dt) * (w.T - Tf)
        d = w.disp_mm(ev.SIG)
        d0 = d if d0 is None else d0
        out.append((Tf, d0 - d, p.resistance(1.0 - w.xi)))      # T1 seen, stroke (mm, contraction positive), R
        ev.plant_step(w, float(I[k]), dt, ev.SIG, T_amb)
    return np.array(out)


def crossing(T, s, frac):
    smax = s.max()
    if smax <= 0:
        return np.nan
    i = np.argmax(s >= frac * smax)
    return float(T[i])


def validate(t, I, V, T1, d, T_amb):
    best = None
    for tau in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):                 # thermocouple lag: the only fitted number
        sim = simulate(t, I, T_amb, tau=tau)
        e = float(np.sqrt(np.mean((sim[:, 0] - T1) ** 2)))
        if best is None or e < best[0]:
            best = (e, tau, sim)
    e, tau, sim = best
    stroke_meas = d[0] - d
    res = {}
    res["rms_T"] = (e, e <= CRIT["rms_T"])
    ts_m, ts_s = crossing(T1, stroke_meas, 0.10), crossing(sim[:, 0], sim[:, 1], 0.10)
    tf_m, tf_s = crossing(T1, stroke_meas, 0.90), crossing(sim[:, 0], sim[:, 1], 0.90)
    res["dT_start"] = (ts_m - ts_s, abs(ts_m - ts_s) <= CRIT["dT_start"])
    res["dT_finish"] = (tf_m - tf_s, abs(tf_m - tf_s) <= CRIT["dT_finish"])
    sm, ss = float(stroke_meas.max()), float(sim[:, 1].max())
    res["stroke_rel"] = ((sm - ss) / ss, abs(sm - ss) / ss <= CRIT["stroke_rel"])
    on = I > 0.5 * I.max()
    Rm = V[on] / I[on]
    cold = float(np.median(Rm[:max(5, len(Rm) // 20)])), float(np.median(sim[on, 2][:max(5, len(Rm) // 20)]))
    hot = float(np.median(Rm[-max(5, len(Rm) // 20):])), float(np.median(sim[on, 2][-max(5, len(Rm) // 20):]))
    rel = max(abs(cold[0] / cold[1] - 1), abs(hot[0] / hot[1] - 1))
    res["R_rel"] = (rel, rel <= CRIT["R_rel"])
    return res, tau


def report(res, tau):
    print(f"thermocouple lag fitted: {tau} s")
    names = dict(rms_T="RMS temperature error [K]", dT_start="start of stroke, T meas - model [K]", dT_finish="90 % of stroke, T meas - model [K]",
                 stroke_rel="final stroke, relative to model", R_rel="resistance, worst of cold/hot, relative")
    for k, (v, ok) in res.items():
        print(f"  {names[k]:44s} {v:+8.3f}   {'PASS' if ok else 'FAIL'}")
    allok = all(ok for _, ok in res.values())
    print("VALIDATED for this run" if allok else "NOT validated: at least one criterion failed (identify h, m c, R and re-run before touching the controller)")
    return allok


def selftest(perturb):
    """synthetic bench log: heating at constant current, then rest; 'true' wire differs from the model by `perturb` in h and m c"""
    rng = np.random.default_rng(1)
    p = ev.NOM
    pt = replace(p, h=p.h * perturb, c=p.c * perturb)
    t = np.arange(0, 400, 0.05)
    I = np.where((t > 5) & (t < 150), 0.9, 0.0)
    w = ev.make_plant(pt, 22.0, 1.0, ev.SIG)
    T1, V, d, Tf = [], [], [], 22.0
    for k in range(len(t)):
        Tf += 0.05 / (1.0 + 0.05) * (w.T - Tf)
        T1.append(Tf + 0.1 * rng.standard_normal()); V.append(I[k] * p.resistance(1.0 - w.xi) * (1 + 0.005 * rng.standard_normal()))
        d.append(w.disp_mm(ev.SIG) + 0.02 * rng.standard_normal())
        ev.plant_step(w, float(I[k]), 0.05, ev.SIG, 22.0)
    return t, I, np.array(V), np.array(T1), np.array(d)


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        for pert in (1.0, 1.6):
            print(f"== selftest: true h and c x{pert}")
            t, I, V, T1, d = selftest(pert)
            res, tau = validate(t, I, V, T1, d, 22.0)
            report(res, tau)
        sys.exit(0)
    amb = float(sys.argv[sys.argv.index("--ambient") + 1]) if "--ambient" in sys.argv else 22.0
    a = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
    res, tau = validate(a["t"], a["I"], a["V"], a["T1"], a["d"], amb)
    sys.exit(0 if report(res, tau) else 1)
