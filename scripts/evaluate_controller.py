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
from sma_control.cycle_controller import CycleController, CycleCfg, cfg_for_wire, HomingRoutine
from sma_control.spring import SMAWireSpring, BiasSpring
from sma_control.elements import SMAWireElements
from sma_control.safety import SafetyMonitor

D_MM = float(os.environ.get("SMA_D_MM", "2.0"))                 # wire diameter [mm]
I_MAX = float(os.environ["SMA_IMAX"]) if "SMA_IMAX" in os.environ else None       # current limit [A] (default: SMAParams default 35 A)
NOM = SMAParams(d=D_MM * 1e-3, **({} if I_MAX is None else dict(Imax=I_MAX)))
SIG = float(os.environ.get("SMA_SIG", "250.0"))                    # MPa, stress at full extension (preload); controller knows it
TS = 0.05
K_SPRING = float(os.environ["SMA_SPRING_K"]) if "SMA_SPRING_K" in os.environ else None      # N/mm, None = constant stress
OVERSHOOT_EXIT = "SMA_NO_OVERSHOOT_EXIT" not in os.environ
GUARD = "SMA_NO_GUARD" not in os.environ                  # independent safety monitor between controller and plant
HOMING = "SMA_NO_HOMING" not in os.environ                # homing cycle before the move sequence (measures end points and windows)


def make_cfg():
    return cfg_for_wire(NOM, Ts=TS, overshoot_exit=OVERSHOOT_EXIT)


def make_plant(p, T0, xi0, sigma_true, struct=None):
    """struct: dict(spread_K, tail_frac, tail_K) -> parallel-elements plant (needs the spring)"""
    if K_SPRING is None:
        return SMAWireCycle(p, T0=T0, xi0=xi0)
    if struct is not None:
        return SMAWireElements(p, BiasSpring(sigma_true, K_SPRING), T0=T0, xi0=xi0, **struct)
    return SMAWireSpring(p, BiasSpring(sigma_true, K_SPRING), T0=T0, xi0=xi0)


def make_ctl(sigma_ctl, cfg, p=None):
    p = p or NOM
    if K_SPRING is None:
        return CycleController(p, sigma_ctl, cfg)
    return CycleController(p, sigma_ctl, cfg, K_spring=BiasSpring(sigma_ctl, K_SPRING).K(p))


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


TRIPS = []        # (label, reason, time) of every safety-monitor trip in this run


class Guard:
    """safety monitor + independent over-temperature channel T2 + fault injection (fault active from t_fault)
    faults: tc_open | tc_stuck (T1 frozen at its value then) | tc_offset (T1 reads 30 K low) | d_stuck | hang | stuck_on | wire_open"""
    def __init__(self, p, Ts, fault=None, t_fault=0.0, seed=2, label="", T_lag=0.0, protect=True, use_T2=True):
        self.mon, self.Ts, self.p, self.T_lag, self.T2f, self.protect = SafetyMonitor(p, Ts), Ts, p, T_lag, None, protect
        self.use_T2 = use_T2
        self._d_frozen, self.Tmax_after = None, 0.0
        self.rows, self.arm_info = [], None            # replay log for the C port: one row per sample (see c/replay.c)
        self.fault, self.tf, self.label = fault, t_fault, label
        self.rng, self.t, self.beat, self.I_prev, self.last_cmd, self._frozen, self._armed = np.random.default_rng(seed), 0.0, 0, 0.0, 0.0, None, False

    def _on(self):
        return self.fault is not None and self.t >= self.tf

    def sense(self, Tm, d, w):
        """apply sensor faults; returns (T1, d)"""
        if self._on():
            if self.fault == "tc_open":
                Tm = 0.0
            elif self.fault == "tc_stuck":
                self._frozen = Tm if self._frozen is None else self._frozen
                Tm = self._frozen
            elif self.fault == "tc_offset":
                Tm = Tm - 30.0
            elif self.fault == "d_stuck":
                self._d_frozen = d if self._d_frozen is None else self._d_frozen
                d = self._d_frozen
        return Tm, d

    def apply(self, w, I_cmd, T1, d, phase=1, tgt=0.0):
        mon, p = self.mon, self.p
        I_ctl = I_cmd
        if not self._armed:
            T2a = (w.T + 0.2 * self.rng.standard_normal()) if self.use_T2 else None
            mon.arm(T1, T2a)
            self.arm_info = (T1, T2a)
            self._armed = True
        if self._on() and self.fault == "hang":
            I_cmd = self.last_cmd                                  # controller frozen: command and heartbeat stop
        else:
            self.last_cmd, self.beat = I_cmd, self.beat + 1
        self.T2f = w.T if self.T2f is None else self.T2f
        if self.T_lag > 0:
            self.T2f += self.Ts / self.T_lag * (w.T - self.T2f)                # second channel is the same sensor type: same lag
        else:
            self.T2f = w.T
        T2 = self.T2f + 0.2 * self.rng.standard_normal() if self.use_T2 else None
        R = p.resistance(1.0 - w.xi) * (3.0 if (self._on() and self.fault == "wire_open") else 1.0)
        V = self.I_prev * R * (1.0 + 0.01 * self.rng.standard_normal())
        if self.protect:
            Ip, relay = mon.step(I_cmd, self.I_prev, T1, T2, d, V, self.beat)
        else:                                                      # same fault injection, no monitor (baseline)
            Ip, relay = float(np.clip(I_cmd, 0.0, p.Imax)), True
        I = (p.Imax if (self._on() and self.fault == "stuck_on") else Ip) if relay else 0.0
        self.rows.append((phase, tgt, T1, d, I_ctl, I_cmd, self.I_prev, T2, V, self.beat, Ip, relay))
        self.I_prev, self.t = I, self.t + self.Ts
        if self._on():
            self.Tmax_after = max(self.Tmax_after, w.T)
        if self.protect and mon.tripped and not any(x[0] == self.label and x[1] == mon.reason for x in TRIPS):
            TRIPS.append((self.label, mon.reason, self.t))
        return I


def new_guard(cfg, sens=None, **kw):
    return Guard(NOM, cfg.Ts, T_lag=getattr(sens, "T_lag", 0.0), **kw) if GUARD else None


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


def do_homing(w, ctl, sens, cfg, T_amb, sigma_true, t_max=6000.0, guard=None):
    """run the homing cycle on the plant; returns dict(t, Tmax, E, fault, cal, guard)"""
    h, t, Tmax, E0 = HomingRoutine(ctl), 0.0, w.T, w.E_joule
    while not h.done and t < t_max:
        Tm, d = sens.read(w.T, w.disp_mm(sigma_true), cfg.Ts)
        if guard:
            Tm, d = guard.sense(Tm, d, w)
        I = h.step(Tm, d)
        if guard:
            I = guard.apply(w, I, Tm, d, phase=0)
            if guard.protect and guard.mon.tripped:
                h.fault, h.done = "safety_trip:" + guard.mon.reason, True
        plant_step(w, I, cfg.Ts, sigma_true, T_amb)
        t += cfg.Ts
        Tmax = max(Tmax, w.T)
    if h.cal is not None:
        ctl.apply_calibration(h.cal)
        if guard:
            guard.mon.set_calibration((h.cal.eps_ext - h.cal.eps_con) * NOM.L * 1e3, h.cal.Af)
    elif not h.done:
        h.fault = "timeout"
    return dict(t=t, Tmax=Tmax, E=w.E_joule - E0, fault=h.fault, cal=h.cal, guard=guard)


def sequence(targets=(1.0, 0.5, 0.8, 0.0), hold_s=30.0, t_max=5000.0, p_true=None, sigma_true=SIG, sigma_ctl=SIG,
             sens=None, T_amb=22.0, struct=None, guard=None, label="", p_ctl=None):
    """event-driven target sequence: next target after the true r has been within 2 % for hold_s seconds"""
    cfg = make_cfg()
    w = make_plant(p_true or NOM, T_amb, 1.0, sigma_true, struct)
    ctl = make_ctl(sigma_ctl, cfg, p_ctl)
    sens = sens or Sensors(T0=T_amb)
    guard = guard or new_guard(cfg, sens=sens, label=label)
    home = do_homing(w, ctl, sens, cfg, T_amb, sigma_true, guard=guard) if HOMING else dict(t=0.0, Tmax=T_amb, E=0.0, fault=None, cal=None, guard=guard)
    if home["fault"]:
        return np.zeros((1, 7)), [], home
    E0 = w.E_joule
    log, seg, k_t, t_in = [], [], 0, None
    for k in range(int(t_max / cfg.Ts)):
        t = k * cfg.Ts
        tgt = targets[k_t]
        Tm, d = sens.read(w.T, w.disp_mm(sigma_true), cfg.Ts)
        r_true = 1.0 - w.xi
        if guard:
            Tm, d = guard.sense(Tm, d, w)
        I = ctl.step(tgt, Tm, d)
        if guard and cfg.heat_only and guard.mon.stroke_mm is None:           # no homing: the stall check works from the model stroke, then from the learned one
            guard.mon.set_calibration(abs(ctl.st.strain(0.0) - ctl.st.strain(1.0)) * ctl.p.L * 1e3, ctl.T_limit_model)
        if guard and cfg.heat_only and ctl.cal is not None and not getattr(guard, "_cal_done", False):
            guard._cal_done = True
            guard.mon.set_calibration(abs(ctl.cal.eps_ext - ctl.cal.eps_con) * ctl.p.L * 1e3, ctl.cal.Af)
        if guard:
            I = guard.apply(w, I, Tm, d, phase=1, tgt=tgt)
        plant_step(w, I, cfg.Ts, sigma_true, T_amb)
        log.append((t, w.T, 1.0 - w.xi, tgt, I, Tm, w.E_joule - E0))
        if guard and guard.protect and guard.mon.tripped:
            home["trip"] = (guard.mon.reason, t)
            break
        if abs(r_true - tgt) < 0.02 and ctl.mode == "hold":
            t_in = t if t_in is None else t_in
            if t - t_in >= hold_s:
                seg.append((tgt, t_in, t)); k_t += 1; t_in = None
                if k_t == len(targets):
                    break
        elif abs(r_true - tgt) >= 0.02:
            t_in = None
    return np.array(log), seg, home


def seq_summary(L, seg, home=None, n_targets=4):
    home = home or dict(t=0.0, Tmax=0.0, fault=None)
    errs, t_prev = [], 0.0
    for tgt, tr, te in seg:
        m = (L[:, 0] >= tr) & (L[:, 0] <= te)
        errs.append(abs(L[m, 2].mean() - tgt))
    return dict(reached=len(seg), n=n_targets, max_err=100 * max(errs) if errs else np.nan, Tmax=max(L[:, 1].max(), home["Tmax"]),
                t_total=L[-1, 0] if len(seg) < n_targets else seg[-1][2], E=L[-1, 6], home_t=home["t"], home_fault=home["fault"])


def home_note(s):
    return f"   homing FAULT: {s['home_fault']}" if s["home_fault"] else ""


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
            ("ambient 15 C", dict(T_amb=15.0)), ("ambient 26 C", dict(T_amb=26.0)),
            ("ambient 26 C, preload 290", dict(T_amb=26.0, sigma_true=290.0, sigma_ctl=290.0))]


def sec_robust(fig=None):
    print("== PARAMETER ERRORS: sequence r = 1.0 → 0.5 → 0.8 → 0.0 (controller built from nominal numbers) ==")
    print(f"{'case':20s} {'targets reached':>15s} {'max hold err %':>14s} {'Tmax':>6s} {'time [s]':>9s} {'E [J]':>7s} {'homing [s]':>10s}")
    rows = []
    for name, kw in case_list():
        L, seg, home = sequence(**kw, label=name)
        s = seq_summary(L, seg, home)
        rows.append((name, s))
        print(f"{name:20s} {s['reached']:>9d} of {s['n']} {s['max_err']:14.2f} {s['Tmax']:6.1f} {s['t_total']:9.0f} {s['E']:7.0f} {s['home_t']:10.0f}{home_note(s)}")
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


STRUCT_CASES = [("single element (controller's own model)", None),
                ("spread 1 K (window about 4 K wider)", dict(spread_K=1.0)),
                ("spread 2 K (about 7 K wider)", dict(spread_K=2.0)),
                ("spread 3 K (about 10 K wider)", dict(spread_K=3.0)),
                ("spread 5 K (about 16 K wider)", dict(spread_K=5.0)),
                ("spread 2 K + 25 % tail to 20 K", dict(spread_K=2.0, tail_frac=0.25, tail_K=20.0)),
                ("spread 3 K + 25 % tail to 20 K", dict(spread_K=3.0, tail_frac=0.25, tail_K=20.0))]


def sec_struct(fig=None):
    if K_SPRING is None:
        print("structural cases need the spring (SMA_SPRING_K)"); return
    print("== STRUCTURAL MISMATCH: plant = parallel elements with spread of transformation temperatures; controller unchanged ==")
    print(f"{'case':42s} {'targets reached':>15s} {'max hold err %':>14s} {'Tmax':>6s} {'time [s]':>9s} {'E [J]':>7s} {'homing [s]':>10s}")
    rows = []
    for name, st in STRUCT_CASES:
        L, seg, home = sequence(struct=st, label=name)
        s_ = seq_summary(L, seg, home)
        rows.append((name, L, seg))
        print(f"{name:42s} {s_['reached']:>9d} of {s_['n']} {s_['max_err']:14.2f} {s_['Tmax']:6.1f} {s_['t_total']:9.0f} {s_['E']:7.0f} {s_['home_t']:10.0f}{home_note(s_)}")
    if fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        sty(plt)
        f, ax = plt.subplots(2, 4, figsize=(17, 7.5), sharey=True)
        for a, (name, L, seg) in zip(ax.ravel(), rows):
            a.plot(L[:, 0], L[:, 3], color=OR, lw=1.3, ls="--"); a.plot(L[:, 0], L[:, 2], color=INK, lw=1.8)
            a.set_title(name, loc="left", fontsize=8.5); a.set_xlabel("time [s]")
        ax[0, 0].set_ylabel("r"); ax[1, 0].set_ylabel("r")
        f.tight_layout(); f.savefig(Path(fig) / "fig17_structural.png", dpi=150, facecolor=BG)


def sec_sensors(fig=None):
    print("== SENSORS: sequence r = 1.0 → 0.5 → 0.8 → 0.0 ==")
    print(f"{'case':36s} {'targets reached':>15s} {'max hold err %':>14s} {'Tmax':>6s} {'time [s]':>9s} {'E [J]':>7s} {'homing [s]':>10s}")
    for name, kw in (("nominal (0.1 K, 0.02 mm)", {}), ("T noise 0.5 K", dict(Tn=0.5)), ("displacement noise 0.1 mm", dict(dn=0.1)),
                     ("displacement quantised 0.05 mm", dict(quant=0.05)), ("thermocouple lag 1 s", dict(T_lag=1.0)),
                     ("thermocouple lag 3 s", dict(T_lag=3.0)),
                     ("all: 0.3 K, 0.05 mm, q 0.02, lag 1 s", dict(Tn=0.3, dn=0.05, quant=0.02, T_lag=1.0))):
        L, seg, home = sequence(sens=Sensors(**kw), label=name)
        s = seq_summary(L, seg, home)
        print(f"{name:36s} {s['reached']:>9d} of {s['n']} {s['max_err']:14.2f} {s['Tmax']:6.1f} {s['t_total']:9.0f} {s['E']:7.0f} {s['home_t']:10.0f}{home_note(s)}")


def sec_safety(fig=None):
    print("== SAFETY: injected faults, with and without the independent monitor ==")
    cfg = make_cfg()
    print(f"controller over-temperature limit {make_ctl(SIG, cfg).T_limit:.1f} C, absolute {cfg.T_abs_max:.0f} C; monitor hard trip {SafetyMonitor(NOM, cfg.Ts).c.T_trip:.0f} C")
    faults = [("thermocouple open (reads 0 C)", "tc_open", (20.0, 290.0)), ("thermocouple frozen", "tc_stuck", (20.0, 290.0)),
              ("thermocouple reads 30 K low", "tc_offset", (20.0, 290.0)),
              ("position sensor frozen", "d_stuck", (20.0, 250.0)), ("controller hangs, current frozen", "hang", (20.0, 250.0)),
              ("power stage stuck on", "stuck_on", (20.0, 290.0)), ("wire contact fails (R x3)", "wire_open", (20.0, 290.0))]
    T_trip = SafetyMonitor(NOM, cfg.Ts).c.T_trip
    print(f"{'fault (second channel T2 present)':36s} {'at [s]':>6s} | {'no monitor: Tmax after':>22s} | {'stopped by':>28s} {'t [s]':>7s} {'Tmax after':>10s}  result")
    n_ok = n = 0
    rows = [(nm, f, tf, True) for nm, f, tfs in faults for tf in tfs] + [(nm + " [no T2]", f, tf, False) for nm, f, tfs in faults[:3] for tf in tfs]
    for name, f, t_f, use_T2 in rows:
        if not use_T2 and name.endswith("[no T2]") and t_f == 20.0 and f == faults[0][1]:
            print(f"{'fault (single thermocouple, envelope check only)':36s}")
        res = []
        for protect in (False, True):
            g = Guard(NOM, cfg.Ts, fault=f, t_fault=t_f, label=name, protect=protect, use_T2=use_T2)
            L, seg, home = sequence(t_max=900.0, guard=g, sens=Sensors(T0=22.0))
            by = g.mon.reason if (protect and g.mon.tripped) else (("homing: " + home["fault"].replace("safety_trip:", "")) if (protect and home["fault"]) else None)
            res.append((g.Tmax_after, by, g.t))
        (Tn, _, _), (Tp, by, tt) = res
        ok = by is not None and Tp < T_trip + 10.0
        n, n_ok = n + 1, n_ok + ok
        print(f"{name:36s} {t_f:6.0f} | {Tn:22.1f} | {str(by):>28s} {tt:7.1f} {Tp:10.1f}  {'ok' if ok else 'NOT STOPPED'}")
    print(f"{n_ok} of {n} fault runs stopped below the trip temperature + 10 K; fail-safe state = power off, wire extends under the bias spring")
    # no-fault runs with the monitor are in every other section ('safety monitor: no trips in any run')


# ====================================================================== style
BG, INK, MUT, GR, OR, BL, GRN = "#fcfcfb", "#0b0b0b", "#52514e", "#a8a7a2", "#eb6834", "#2a78d6", "#1baf7a"


def sty(plt):
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": "#e6e5e1", "axes.labelcolor": MUT,
                         "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": "#e6e5e1",
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5})


SECTIONS = dict(struct=sec_struct, steps=sec_steps, ramps=sec_ramps, disturb=sec_disturb, robust=sec_robust, sensors=sec_sensors, safety=sec_safety)

if __name__ == "__main__":
    SECTIONS[sys.argv[1]](sys.argv[2] if len(sys.argv) > 2 else None)
    if GUARD and sys.argv[1] != "safety":
        print("safety monitor: " + ("no trips in any run" if not TRIPS else "TRIPS " + "; ".join(f"{a}: {b} at {c:.0f} s" for a, b, c in TRIPS)))
