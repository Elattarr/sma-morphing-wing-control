"""Demo: PID (+ static feedforward) on the 2 mm thesis wire, plus a 0.5 mm comparison.
Run:  python run_demo.py      -> writes sma_demo_*.png and prints metrics"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sma_wire import SMAParams, SMAWire, PIDGains, simulate, step_metrics, make_ff

# ---------------------------------------------------------------- 2 mm wire
p = SMAParams(); p.Imax = 10.0
ref = lambda t: 0 if t < 5 else (4.0 if t < 150 else (6.0 if t < 300 else 3.0))
ff = make_ff(p)
cases = {
    "PID":        (PIDGains(Kp=9.7, Ki=1.7, Kd=0, Kt=1.0), None),
    "PID + FF":   (PIDGains(Kp=23.7, Ki=7.2, Kd=0, Kt=1.0), ff),
}
fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
for name, (g, f) in cases.items():
    o = simulate(p, g, ref, 450, dt=0.02, ff_fn=f)
    ax[0].plot(o["t"], o["ytrue"], label=name)
    ax[1].plot(o["t"], o["I"], label=name)
    ax[2].plot(o["t"], o["T"], label=name)
    print(name)
    for ts, rf, y0, te in [(5, 4, 0, 150), (150, 6, 4, 300), (300, 3, 6, 450)]:
        m = (o["t"] >= ts) & (o["t"] < te)
        mt = step_metrics(o["t"][m], o["ytrue"][m], rf, ts, y0)
        print(f"  step to {rf} mm:", {k: round(float(v), 1) for k, v in mt.items()},
              "Tmax %.0f C" % o["T"][m].max())
ax[0].plot(o["t"], o["r"], "k--", label="reference")
ax[0].set_ylabel("stroke [mm]"); ax[1].set_ylabel("current [A]"); ax[2].set_ylabel("wire T [C]")
ax[2].axhline(p.Af + p.sigma0 / p.beta, ls=":", c="gray"); ax[2].text(5, p.Af + p.sigma0 / p.beta + 2, "Af at bias stress")
ax[2].set_xlabel("t [s]")
for a in ax: a.grid(alpha=.3); a.legend(loc="best")
ax[0].set_title("2 mm x 200 mm wire, closed loop (heating fast, cooling passive and slow)")
fig.tight_layout(); fig.savefig("sma_demo_steps.png", dpi=130)

# ----------------------------------------------- hysteresis loop (quasi-static)
w = SMAWire(SMAParams(noise_mm=0)); T_, s_ = [], []
for Tt in list(np.linspace(22, 95, 500)) + list(np.linspace(95, 22, 500)):
    w.step(0.0, 0.01, T_force=Tt)
    T_.append(Tt); s_.append(w.stroke_mm())
fig, a = plt.subplots(figsize=(6, 4.2))
a.plot(T_, s_); a.set_xlabel("T [C]"); a.set_ylabel("stroke [mm]"); a.grid(alpha=.3)
a.set_title("Stroke vs temperature at bias stress (stress-shifted hysteresis)")
fig.tight_layout(); fig.savefig("sma_demo_hysteresis.png", dpi=130)

# --------------------------------------- why wire diameter matters for control
fig, a = plt.subplots(figsize=(7, 4))
for d, Imax in [(2.0e-3, 10.0), (0.5e-3, 3.0)]:
    q = SMAParams(d=d, Imax=Imax, noise_mm=0)
    tau = q.mass * q.c / (q.h * q.As_surf)
    w = SMAWire(q); t, T = [], []
    for k in range(int(160 / 0.02)):
        # heat at 1.6 x hold-current estimate, then 0
        I = Imax if k * 0.02 < 15 else 0.0
        w.step(I, 0.02); t.append(k * 0.02); T.append(w.T)
    a.plot(t, T, label=f"d = {d*1e3:.1f} mm  (tau_cool = {tau:.0f} s)")
a.set_xlabel("t [s]"); a.set_ylabel("T [C]"); a.grid(alpha=.3); a.legend()
a.set_title("Heat 15 s then cool: bandwidth set by wire diameter")
fig.tight_layout(); fig.savefig("sma_demo_diameter.png", dpi=130)
