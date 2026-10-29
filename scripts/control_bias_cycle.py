"""Closed-loop position control of the one-way wire under a constant bias stress (both directions).
python scripts/control_bias_cycle.py [figure_dir]

Plant      : SMAWireCycle (MEFF, one-way), sigma_true, ambient 22 C, start xi = 1 at 22 C (stored, extended)
Controller : CycleController with sigma_ctl (what the controller believes)
Sensors    : thermocouple noise 0.1 K, displacement noise 0.02 mm
Targets    : r = 1.0 (full contraction) -> 0.5 -> 0.8 -> 0.0, each held at least 30 s once reached
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sma_control import SMAParams
from sma_control.cycle_plant import SMAWireCycle
from sma_control.cycle_controller import CycleController, CycleCfg
from sma_control.spring import SMAWireSpring, BiasSpring

TARGETS = [1.0, 0.5, 0.8, 0.0]
HOLD_S, T_AMB, TN = 30.0, 22.0, 0.1


def run(sigma_true=250.0, sigma_ctl=250.0, targets=TARGETS, t_max=8000.0, seed=1, k_true=None, k_ctl=None):
    """k_true = None: constant stress sigma_true.  k_true = spring rate [N/mm]: bias spring, sigma_true = stress at full extension;
    k_ctl = spring rate the controller believes (default: the true one)"""
    p = SMAParams()
    cfg = CycleCfg()
    if k_true is None:
        w = SMAWireCycle(p, T0=T_AMB, xi0=1.0)
        ctl = CycleController(p, sigma_ctl, cfg)
    else:
        w = SMAWireSpring(p, BiasSpring(sigma_true, k_true), T0=T_AMB, xi0=1.0)
        ctl = CycleController(p, sigma_ctl, cfg, K_spring=BiasSpring(sigma_ctl, k_true if k_ctl is None else k_ctl).K(p))
    rng = np.random.default_rng(seed)
    log, seg, k_t, t_reached = [], [], 0, None
    n = int(t_max / cfg.Ts)
    for k in range(n):
        t = k * cfg.Ts
        r_true = 1.0 - w.xi
        tgt = targets[k_t]
        d_meas = w.disp_mm(sigma_true) + p.noise_mm * rng.standard_normal()
        I = ctl.step(tgt, w.T + TN * rng.standard_normal(), d_meas)
        E0 = w.E_joule
        w.step(I, cfg.Ts, sigma_true, T_AMB) if k_true is None else w.step(I, cfg.Ts, T_AMB)
        log.append((t, w.T, 1.0 - w.xi, tgt, I, ctl.mode == "hold"))
        if abs(r_true - tgt) < 0.02 and ctl.mode == "hold":
            t_reached = t if t_reached is None else t_reached
            if t - t_reached >= HOLD_S:
                seg.append((tgt, t_reached, t, w.E_joule))
                k_t += 1; t_reached = None
                if k_t == len(targets):
                    break
        elif t_reached is not None and abs(r_true - tgt) >= 0.02:
            t_reached = None
    return np.array(log), seg, ctl


def report(L, seg, ctl, label):
    print(f"-- {label} --  (controller windows: A_s {ctl.As:.1f}, A_f {ctl.Af:.1f}, forward start {ctl.Tfs:.1f}, "
          f"T_safe {ctl.T_safe:.1f} C; max extension rate {ctl.rate_dn:.4f} /s)")
    t_prev, E_prev = 0.0, 0.0
    for tgt, tr, te, E in seg:
        m = (L[:, 0] >= tr) & (L[:, 0] <= te)
        err = L[m, 2].mean() - tgt
        win = (L[:, 0] >= t_prev) & (L[:, 0] <= te)
        r_seg = L[win, 2]
        print(f"target {tgt:.1f}: reached after {tr - t_prev:6.1f} s (|err| < 2 %), mean hold error {err:+.3f}, "
              f"T during hold {L[m, 1].mean():5.1f} C, Tmax in segment {L[win, 1].max():5.1f} C, "
              f"energy {E - E_prev:6.0f} J")
        t_prev, E_prev = te, E
    if len(seg) < len(TARGETS):
        print(f"  !! only {len(seg)} of {len(TARGETS)} targets reached; final r = {L[-1, 2]:.3f}, T = {L[-1, 1]:.1f} C")


def figure(L, outdir, name="fig9_bias_cycle_control.png", title="Bidirectional position control, 250 MPa bias (controller knows sigma)"):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BG, INK, MUT, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": GRID, "axes.labelcolor": MUT,
                         "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9.5})
    fig, ax = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)
    ax[0].plot(L[:, 0], L[:, 3], color="#eb6834", lw=1.6, ls="--", label="target r")
    ax[0].plot(L[:, 0], L[:, 2], color=INK, lw=2, label="r (plant)")
    ax[0].set_ylabel("recovered fraction r"); ax[0].legend(frameon=False, loc="center right")
    ax[1].plot(L[:, 0], L[:, 1], color="#2a78d6", lw=1.6)
    ax[1].set_ylabel("T [°C]")
    ax[2].plot(L[:, 0], L[:, 4], color="#1baf7a", lw=1.2)
    ax[2].set_ylabel("current I [A]"); ax[2].set_xlabel("time [s]")
    ax[0].set_title(title, loc="left", fontsize=10)
    fig.tight_layout(); Path(outdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(outdir) / name, dpi=160, facecolor=BG)


if __name__ == "__main__":
    L, seg, ctl = run()
    report(L, seg, ctl, "nominal: sigma_true = sigma_ctl = 250 MPa")
    for st, sc in ((250.0, 225.0), (250.0, 275.0)):
        L2, seg2, ctl2 = run(st, sc)
        report(L2, seg2, ctl2, f"sigma_true {st:.0f}, controller believes {sc:.0f} MPa")
    print()
    Ls, segs, ctls = run(250.0, 250.0, k_true=20.0)
    report(Ls, segs, ctls, "BIAS SPRING 20 N/mm, 250 MPa preload at full extension (controller knows the spring)")
    for st, sc, kt, kc in ((250.0, 230.0, 20.0, 20.0), (250.0, 270.0, 20.0, 20.0), (250.0, 250.0, 20.0, 15.0), (250.0, 250.0, 20.0, 25.0)):
        L2, seg2, ctl2 = run(st, sc, k_true=kt, k_ctl=kc)
        report(L2, seg2, ctl2, f"spring: preload true {st:.0f} / believed {sc:.0f} MPa, rate true {kt:.0f} / believed {kc:.0f} N/mm")
    if len(sys.argv) > 1:
        figure(L, sys.argv[1])
        figure(Ls, sys.argv[1], "fig16_spring_cycle_control.png", "Bidirectional position control with a bias spring (20 N/mm, 250 MPa at full extension)")
