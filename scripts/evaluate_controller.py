"""Controller evaluation on our model: one-way wire (MEFF), 250 MPa bias, bidirectional position controller.

python scripts/evaluate_controller.py steps|ramps|disturb|robust|sensors|safety [figure_dir]
Environment: SMA_D_MM=<mm> wire diameter (default 2), SMA_IMAX=<A> current limit (default 35)
             SMA_SPRING_K=<N/mm>  bias spring instead of a constant stress (SIG = stress at full extension, controller knows the spring)
             SMA_NO_OVERSHOOT_EXIT=1  switch the overshoot exit of the controller off (comparison)

Plant      : SMAWireCycle (true parameters, possibly different from the controller's)
Controller : CycleController built from the NOMINAL SMAParams and sigma_ctl, Ts = 0.05 s
Metrics    : see each section. r = recovered fraction (0 extended ... 1 contracted); errors are in % of stroke (r x 100).
"""
import os
import sys
from dataclasses import replace
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sma_control import SMAParams
from sma_control.cycle_plant import SMAWireCycle
from sma_control.cycle_controller import CycleController, CycleCfg, cfg_for_wire
from sma_control.spring import SMAWireSpring, BiasSpring

D_MM = float(os.environ.get("SMA_D_MM", "2.0"))                 # wire diameter [mm]
I_MAX = float(os.environ["SMA_IMAX"]) if "SMA_IMAX" in os.environ else None       # current limit [A] (default: SMAParams default 35 A)
NOM = SMAParams(d=D_MM * 1e-3, **({} if I_MAX is None else dict(Imax=I_MAX)))
SIG = 250.0
TS = 0.05
K_SPRING = float(os.environ["SMA_SPRING_K"]) if "SMA_SPRING_K" in os.environ else None      # N/mm, None = constant stress
OVERSHOOT_EXIT = "SMA_NO_OVERSHOOT_EXIT" not in os.environ


def make_cfg():
    return cfg_for_wire(NOM, Ts=TS, overshoot_exit=OVERSHOOT_EXIT)


def make_plant(p, T0, xi0, sigma_true):
    if K_SPRING is None:
        return SMAWireCycle(p, T0=T0, xi0=xi0)
    return SMAWireSpring(p, BiasSpring(sigma_true, K_SPRING), T0=T0, xi0=xi0)


def make_ctl(sigma_ctl, cfg):
    if K_SPRING is None:
        return CycleController(NOM, sigma_ctl, cfg)
    return CycleController(NOM, sigma_ctl, cfg, K_spring=BiasSpring(sigma_ctl, K_SPRING).K(NOM))


def plant_step(w, I, dt, sigma_true, T_inf, h_scale=1.0):
    if K_SPRING is None:
        return w.step(I, dt, sigma_true, T_inf, h_scale)
    return w.step(I, dt, T_inf, h_scale)


T_SAFE = make_ctl(SIG, make_cfg()).T_safe          # hold temperature (31.7 C constant stress, 37.7 C with the 20 N/mm spring)
DOCS = Path(__file__).resolve().parents[1] / "docs"


class Sensors:
    """thermocouple and displacement sensor with optional lag, noise, quantisation, failure"""
    def __init__(self, Tn=0.1, dn=0.02, quant=0.0, T_lag=0.0, T_open=False, d_stuck=None, seed=1, T0=22.0):
        self.Tn, self.dn, self.quant, self.T_lag, self.T_open, self.d_stuck = Tn, dn, quant, T_lag, T_open, d_stuck
        self.rng, self.Tf = np.random.default_rng(seed), T0

    def read(self, T_true, d_true, dt):
        if self.T_lag > 0:
            self.Tf += dt / self.T_lag * (T_true - self.Tf); Tm = self.Tf
        else:
            Tm = T_true
        Tm = 0.0 if self.T_open else Tm + self.Tn * self.rng.standard_normal()
        d = d_true if self.d_stuck is None else self.d_stuck
        d = d + self.dn * self.rng.standard_normal()
        if self.quant > 0:
            d = round(d / self.quant) * self.quant
        return Tm, d


def simulate(r_ref, T_end, p_true=None, sigma_true=SIG, sigma_ctl=SIG, cfg=None, xi0=1.0, T0=22.0,
             T_amb=lambda t: 22.0, h_scale=lambda t: 1.0, sens=None):
    """returns log [t, T, r, r_ref, I, T_meas, E_joule]"""
    cfg = cfg or make_cfg()
    w = make_plant(p_true or NOM, T0, xi0, sigma_true)
    ctl = make_ctl(sigma_ctl, cfg)
    sens = sens or Sensors(T0=T0)
    n = int(T_end / cfg.Ts)
    log = np.zeros((n, 7))
    for k in range(n):
        t = k * cfg.Ts
        Tm, d = sens.read(w.T, w.disp_mm(sigma_true), cfg.Ts)
        I = ctl.step(r_ref(t), Tm, d)
        plant_step(w, I, cfg.Ts, sigma_true, T_amb(t), h_scale(t))
        log[k] = (t, w.T, 1.0 - w.xi, r_ref(t), I, Tm, w.E_joule)
    return log


def sequence(targets=(1.0, 0.5, 0.8, 0.0), hold_s=30.0, t_max=5000.0, p_true=None, sigma_true=SIG, sigma_ctl=SIG,
             sens=None, T_amb=22.0):
    """event-driven target sequence: next target after the true r has been within 2 % for hold_s seconds"""
    cfg = make_cfg()
    w = make_plant(p_true or NOM, T_amb, 1.0, sigma_true)
    ctl = make_ctl(sigma_ctl, cfg)
    sens = sens or Sensors(T0=T_amb)
    log, seg, k_t, t_in = [], [], 0, None
    for k in range(int(t_max / cfg.Ts)):
        t = k * cfg.Ts
        tgt = targets[k_t]
        Tm, d = sens.read(w.T, w.disp_mm(sigma_true), cfg.Ts)
        r_true = 1.0 - w.xi
        I = ctl.step(tgt, Tm, d)
        plant_step(w, I, cfg.Ts, sigma_true, T_amb)
        log.append((t, w.T, 1.0 - w.xi, tgt, I, Tm, w.E_joule))
        if abs(r_true - tgt) < 0.02 and ctl.mode == "hold":
            t_in = t if t_in is None else t_in
            if t - t_in >= hold_s:
                seg.append((tgt, t_in, t)); k_t += 1; t_in = None
                if k_t == len(targets):
                    break
        elif abs(r_true - tgt) >= 0.02:
            t_in = None
    return np.array(log), seg


def seq_summary(L, seg, n_targets=4):
    errs, t_prev = [], 0.0
    for tgt, tr, te in seg:
        m = (L[:, 0] >= tr) & (L[:, 0] <= te)
        errs.append(abs(L[m, 2].mean() - tgt))
    return dict(reached=len(seg), n=n_targets, max_err=100 * max(errs) if errs else np.nan, Tmax=L[:, 1].max(),
                t_total=L[-1, 0] if len(seg) < n_targets else seg[-1][2], E=L[-1, 6])


# ====================================================================== metrics
def step_metrics(L, r0, r1, tail=20.0):
    t, r = L[:, 0], L[:, 2]
    s = r1 - r0
    prog = (r - r0) / s
    def first(level):
        i = int(np.argmax(prog >= level))
        return t[i] if prog[i] >= level else np.nan
    out = np.where(np.abs(r - r1) > 0.02)[0]
    ts = np.nan if (len(out) and out[-1] + 1 >= len(t)) else (t[out[-1] + 1] if len(out) else t[0])
    i_s = int(np.searchsorted(t, ts)) if not np.isnan(ts) else len(t) - 1
    return dict(delay=first(0.02), rise=first(0.9) - first(0.1), settle=ts, os=100 * max(prog.max() - 1.0, 0.0),
                sse=100 * (r[t >= t[-1] - tail].mean() - r1), Tmax=L[:, 1].max(), Imax=L[:, 4].max(),
                t35=(L[:, 4] > 34.9).sum() * TS, E=L[min(i_s, len(t) - 1), 6], Istd=L[t >= t[-1] - tail, 4].std())


def tracking_metrics(L, t_skip=0.0):
    t, r, rr = L[:, 0], L[:, 2], L[:, 3]
    m = t >= t_skip
    e = r[m] - rr[m]
    lags = np.arange(0, 700, 10.0)
    rms = []
    for lag in lags:
        j = np.searchsorted(t, t - lag)
        rms.append(np.sqrt(np.mean((r[m] - rr[np.clip(j, 0, len(t) - 1)][m]) ** 2)))
    return dict(rms=100 * np.sqrt(np.mean(e ** 2)), emax=100 * np.abs(e).max(), lag=lags[int(np.argmin(rms))], rms_lag=100 * min(rms))


def time_to(L, r1, level=0.02):
    i = np.where(np.abs(L[:, 2] - r1) < level)[0]
    return L[i[0], 0] if len(i) else np.nan


# ====================================================================== sections
def sec_steps(fig=None):
    out, runs = [], {}
    print("== STEP RESPONSES (250 MPa bias, noise 0.1 K / 0.02 mm) ==")
    print(f"{'step':>12s} {'delay':>7s} {'rise':>7s} {'settle':>7s} {'OS %':>6s} {'sse %':>7s} {'Tmax':>6s} {'Imax':>5s} {'t@35A':>6s} {'E[J]':>6s} {'I std':>6s}")
    for kind, r0, targets, xi0, T0, T_end in (("contract", 0.0, (0.25, 0.5, 0.75, 1.0), 1.0, 22.0, 120.0),
                                               ("extend", 1.0, (0.75, 0.5, 0.25, 0.0), 0.0, T_SAFE, 1000.0),
                                               ("reversal hot", 1.0, (0.5,), 0.0, 88.0, 1200.0)):
        for r1 in targets:
            L = simulate(lambda t, r1=r1: r1, T_end, xi0=xi0, T0=T0)
            m = step_metrics(L, r0, r1)
            runs[(kind, r1)] = L
            print(f"{kind[:8]:>8s} {r0:.2f}→{r1:.2f} {m['delay']:7.1f} {m['rise']:7.1f} {m['settle']:7.1f} {m['os']:6.2f} {m['sse']:+7.2f} "
                  f"{m['Tmax']:6.1f} {m['Imax']:5.1f} {m['t35']:6.1f} {m['E']:6.0f} {m['Istd']:6.3f}")
    if fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        sty(plt); fig_dir = Path(fig)
        f, ax = plt.subplots(2, 2, figsize=(13, 8.5))
        for r1 in (0.25, 0.5, 0.75, 1.0):
            L = runs[("contract", r1)]; ax[0, 0].plot(L[:, 0], L[:, 2], lw=2, label=f"target {r1}"); ax[0, 0].axhline(r1, color=GR, lw=1, ls="--")
            ax[1, 0].plot(L[:, 0], L[:, 4], lw=1.4)
        ax[0, 0].set_xlim(0, 25); ax[1, 0].set_xlim(0, 25)
        ax[0, 0].set_ylabel("recovered fraction r"); ax[0, 0].set_title("Contraction steps from r = 0 (dashed: targets)", loc="left", fontsize=10); ax[0, 0].legend(frameon=False, fontsize=8.5)
        ax[1, 0].set_ylabel("current I [A]"); ax[1, 0].set_xlabel("time [s]")
        for r1 in (0.75, 0.5, 0.25, 0.0):
            L = runs[("extend", r1)]; ax[0, 1].plot(L[:, 0], L[:, 2], lw=2, label=f"target {r1}"); ax[0, 1].axhline(r1, color=GR, lw=1, ls="--")
            ax[1, 1].plot(L[:, 0], L[:, 1], lw=1.4)
        ax[0, 1].set_ylabel("recovered fraction r"); ax[0, 1].set_title(f"Extension steps from r = 1 (start at {T_SAFE:.1f} °C)", loc="left", fontsize=10); ax[0, 1].legend(frameon=False, fontsize=8.5)
        ax[1, 1].set_ylabel("temperature T [°C]"); ax[1, 1].set_xlabel("time [s]")
        f.tight_layout(); f.savefig(fig_dir / "fig12_step_responses.png", dpi=150, facecolor=BG)


def sec_ramps(fig=None):
    P1, P2 = 2400.0, 900.0
    cases = [("up ramp 0→1 in 300 s", lambda t: float(np.clip(t / 300.0, 0, 1)), 450.0, 1.0, 22.0, 0.0),
             ("down ramp 1→0 in 900 s", lambda t: float(np.clip(1 - t / 900.0, 0, 1)), 1100.0, 0.0, T_SAFE, 0.0),
             (f"sine 0.5±0.2, period {P1:.0f} s", lambda t: 0.5 + 0.2 * np.sin(2 * np.pi * t / P1), 2 * P1, 0.5, T_SAFE, 100.0),
             (f"sine 0.5±0.2, period {P2:.0f} s", lambda t: 0.5 + 0.2 * np.sin(2 * np.pi * t / P2), 3 * P2, 0.5, T_SAFE, 100.0)]
    print("== RAMP AND SINE TRACKING ==")
    print(f"{'case':34s} {'RMS err %':>9s} {'max err %':>9s} {'best lag s':>10s} {'RMS after lag %':>16s} {'Tmax':>6s} {'E[J]':>6s}")
    logs = []
    for name, rf, T_end, xi_r, T0, skip in cases:
        L = simulate(rf, T_end, xi0=xi_r, T0=T0)
        m = tracking_metrics(L, t_skip=skip)
        logs.append((name, L, m))
        print(f"{name:34s} {m['rms']:9.2f} {m['emax']:9.2f} {m['lag']:10.0f} {m['rms_lag']:16.2f} {L[:, 1].max():6.1f} {L[-1, 6]:6.0f}")
    if fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        sty(plt)
        f, ax = plt.subplots(2, 2, figsize=(13, 8.5))
        for a, (name, L, m) in zip(ax.ravel(), logs):
            a.plot(L[:, 0], L[:, 3], color=OR, lw=1.6, ls="--", label="reference"); a.plot(L[:, 0], L[:, 2], color=INK, lw=2, label="r (plant)")
            a.set_title(f"{name}: RMS {m['rms']:.1f} %, lag {m['lag']:.0f} s", loc="left", fontsize=9.5); a.set_xlabel("time [s]"); a.set_ylabel("r")
            a.legend(frameon=False, fontsize=8.5)
        f.tight_layout(); f.savefig(Path(fig) / "fig13_ramp_sine_tracking.png", dpi=150, facecolor=BG)


def sec_disturb(fig=None):
    print("== DISTURBANCES DURING A MOVE ==")
    T_ON = 150.0 * NOM.cooling_tau / SMAParams().cooling_tau        # disturbance starts at about 1/3 of the baseline extension (150 s for the 2 mm wire)
    base_ext = simulate(lambda t: 0.5, 900.0, xi0=0.0, T0=T_SAFE)
    base_con = simulate(lambda t: 1.0, 60.0, xi0=1.0, T0=22.0)
    print(f"disturbance onset: {T_ON:.0f} s")
    print(f"baseline extension 1→0.5: reaches 0.5 (|e|<2 %) at {time_to(base_ext, 0.5):.0f} s, final error {100 * (base_ext[-1, 2] - 0.5):+.2f} %")
    print(f"baseline contraction 0→1: reaches 1 at {time_to(base_con, 1.0):.1f} s")
    cases = [("ext: air flow x2 from onset", dict(h_scale=lambda t: 2.0 if t > T_ON else 1.0), base_ext, 0.5, 900.0, 0.0, T_SAFE),
             ("ext: air flow x0.5 from onset", dict(h_scale=lambda t: 0.5 if t > T_ON else 1.0), base_ext, 0.5, 900.0, 0.0, T_SAFE),
             ("ext: ambient +4 K from onset", dict(T_amb=lambda t: 26.0 if t > T_ON else 22.0), base_ext, 0.5, 900.0, 0.0, T_SAFE),
             ("ext: ambient +8 K from onset", dict(T_amb=lambda t: 30.0 if t > T_ON else 22.0), base_ext, 0.5, 900.0, 0.0, T_SAFE),
             ("con: air flow x3 from 0 s", dict(h_scale=lambda t: 3.0), base_con, 1.0, 60.0, 1.0, 22.0),
             ("con: ambient -10 K", dict(T_amb=lambda t: 12.0), base_con, 1.0, 60.0, 1.0, 12.0)]
    print(f"{'case':34s} {'time to target':>14s} {'final err %':>11s} {'Tmax':>6s} {'E[J]':>6s}   note")
    runs = []
    for name, kw, base, r1, T_end, xi0, T0 in cases:
        L = simulate(lambda t, r1=r1: r1, T_end, xi0=xi0, T0=T0, **kw)
        tt = time_to(L, r1)
        note = "target never reached" if np.isnan(tt) else ""
        print(f"{name:34s} {tt:14.1f} {100 * (L[-1, 2] - r1):+11.2f} {L[:, 1].max():6.1f} {L[-1, 6]:6.0f}   {note}")
        runs.append((name, L, base, r1))
    if fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        sty(plt)
        f, ax = plt.subplots(2, 3, figsize=(15, 8))
        for a, (name, L, base, r1) in zip(ax.ravel(), runs):
            a.plot(base[:, 0], base[:, 2], color=GR, lw=1.6, ls="--", label="no disturbance")
            a.plot(L[:, 0], L[:, 2], color=INK, lw=2, label="disturbed"); a.axhline(r1, color=OR, lw=1, ls=":")
            a.set_title(name, loc="left", fontsize=9.5); a.set_xlabel("time [s]"); a.set_ylabel("r"); a.legend(frameon=False, fontsize=8.5)
            if "con" in name: a.set_xlim(0, 40)
        f.tight_layout(); f.savefig(Path(fig) / "fig14_disturbances.png", dpi=150, facecolor=BG)


def case_list():
    lam = 0.0727                                         # J/mm3, Clausius-Clapeyron value (46 J)
    m = NOM.meff
    wide = replace(m, h=12.0 * m.beta / m.eps_L)
    return [("nominal", {}), ("Lambda = 46 J", dict(p_true=replace(NOM, latent_J_per_mm3=lam))),
            ("c = 320", dict(p_true=replace(NOM, c=320.0))), ("h_c = 5", dict(p_true=replace(NOM, h=5.0))),
            ("h_c = 82", dict(p_true=replace(NOM, h=82.0))),
            ("rho reversed", dict(p_true=replace(NOM, rho_eM=1.0e-6, rho_eA=0.87e-6))),
            ("rho +20 %", dict(p_true=replace(NOM, rho_eM=0.96e-6, rho_eA=1.2e-6))),
            ("rho -20 %", dict(p_true=replace(NOM, rho_eM=0.64e-6, rho_eA=0.8e-6))),
            ("window 12 K", dict(p_true=replace(NOM, meff=wide))),
            ("sigma_ctl -10 %", dict(sigma_ctl=225.0)), ("sigma_ctl +10 %", dict(sigma_ctl=275.0)),
            ("sigma_ctl -20 %", dict(sigma_ctl=200.0)), ("sigma_ctl +20 %", dict(sigma_ctl=300.0)),
            ("ambient 15 C", dict(T_amb=15.0)), ("ambient 26 C", dict(T_amb=26.0))]


def sec_robust(fig=None):
    print("== PARAMETER ERRORS: sequence r = 1.0 → 0.5 → 0.8 → 0.0 (controller built from nominal numbers) ==")
    print(f"{'case':20s} {'targets reached':>15s} {'max hold err %':>14s} {'Tmax':>6s} {'time [s]':>9s} {'E [J]':>7s}")
    rows = []
    for name, kw in case_list():
        L, seg = sequence(**kw)
        s = seq_summary(L, seg)
        rows.append((name, s))
        print(f"{name:20s} {s['reached']:>9d} of {s['n']} {s['max_err']:14.2f} {s['Tmax']:6.1f} {s['t_total']:9.0f} {s['E']:7.0f}")
    if fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        sty(plt)
        f, ax = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
        names = [r[0] for r in rows]; y = np.arange(len(rows))[::-1]
        ax[0].barh(y, [r[1]["max_err"] for r in rows], color=BL); ax[0].set_xlabel("max hold error [% of stroke]")
        ax[1].barh(y, [r[1]["t_total"] for r in rows], color=OR); ax[1].set_xlabel("time for the 4-move sequence [s]")
        ax[2].barh(y, [r[1]["Tmax"] for r in rows], color=GRN); ax[2].set_xlabel("peak temperature [°C]")
        ax[0].set_yticks(y); ax[0].set_yticklabels(names)
        f.tight_layout(); f.savefig(Path(fig) / "fig15_robustness.png", dpi=150, facecolor=BG)


def sec_sensors(fig=None):
    print("== SENSORS: sequence r = 1.0 → 0.5 → 0.8 → 0.0 ==")
    print(f"{'case':36s} {'targets reached':>15s} {'max hold err %':>14s} {'Tmax':>6s} {'time [s]':>9s} {'E [J]':>7s}")
    for name, kw in (("nominal (0.1 K, 0.02 mm)", {}), ("T noise 0.5 K", dict(Tn=0.5)), ("displacement noise 0.1 mm", dict(dn=0.1)),
                     ("displacement quantised 0.05 mm", dict(quant=0.05)), ("thermocouple lag 1 s", dict(T_lag=1.0)),
                     ("thermocouple lag 3 s", dict(T_lag=3.0)),
                     ("all: 0.3 K, 0.05 mm, q 0.02, lag 1 s", dict(Tn=0.3, dn=0.05, quant=0.02, T_lag=1.0))):
        L, seg = sequence(sens=Sensors(**kw))
        s = seq_summary(L, seg)
        print(f"{name:36s} {s['reached']:>9d} of {s['n']} {s['max_err']:14.2f} {s['Tmax']:6.1f} {s['t_total']:9.0f} {s['E']:7.0f}")


def sec_safety(fig=None):
    print("== SAFETY ==")
    A_f = make_ctl(SIG, make_cfg()).Af
    print(f"true A_f (end of recovery, highest stress): {A_f:.1f} C; controller over-temperature limit: {make_ctl(SIG, make_cfg()).T_limit:.1f} C")
    cases = [("thermocouple open (reads 0 C)", dict(sens=Sensors(T_open=True)), SIG),
             ("displacement sensor stuck at start", dict(sens=Sensors(d_stuck=NOM.L * 1e3 * np.sqrt(2 / 3) * NOM.meff.eps_L + SIG / 46000.0 * 200.0)), SIG),
             ("true stress 300, controller assumes 250", dict(sigma_true=300.0), SIG),
             ("true stress 350, controller assumes 250", dict(sigma_true=350.0), SIG)]
    print(f"{'case':42s} {'Tmax true':>10s} {'final r':>8s} {'I at end':>9s}  note")
    for name, kw, _ in cases:
        L = simulate(lambda t: 1.0, 200.0, **kw)
        note = "NO INDEPENDENT PROTECTION: needs a hardware limit" if L[:, 1].max() > 120 else ""
        print(f"{name:42s} {L[:, 1].max():10.1f} {L[-1, 2]:8.3f} {L[-1, 4]:9.1f}  {note}")
    print(f"current limit: {NOM.Imax:.0f} A (clip in plant and controller); repeated cycles: identical by construction (MEFF has no drift)")


# ====================================================================== style
BG, INK, MUT, GR, OR, BL, GRN = "#fcfcfb", "#0b0b0b", "#52514e", "#a8a7a2", "#eb6834", "#2a78d6", "#1baf7a"


def sty(plt):
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": "#e6e5e1", "axes.labelcolor": MUT,
                         "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": "#e6e5e1",
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5})


SECTIONS = dict(steps=sec_steps, ramps=sec_ramps, disturb=sec_disturb, robust=sec_robust, sensors=sec_sensors, safety=sec_safety)

if __name__ == "__main__":
    SECTIONS[sys.argv[1]](sys.argv[2] if len(sys.argv) > 2 else None)
