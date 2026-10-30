"""Log closed-loop Python runs (controller + homing + safety monitor) as replay scenarios for the C port (c/replay.c).

python scripts/export_replay.py <out_dir>      writes <out_dir>/<name>.txt for every scenario in SCENARIOS
Each file: wire / meff / cfg / ctl / safety constants, the arming values of the monitor, then one row per sample:
  phase(0 homing, 1 control)  r_target  T1  d  I_ctl(py)  I_to_monitor  I_meas  have_T2 T2  have_V V  have_d  beat  I_permitted(py)  relay(py)
"""
import os
import sys
from dataclasses import asdict, replace
from pathlib import Path
import numpy as np
os.environ.setdefault("SMA_D_MM", "0.5"); os.environ.setdefault("SMA_IMAX", "3"); os.environ.setdefault("SMA_SPRING_K", "1.25")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_controller as ev
from sma_control.loading import UniaxialMEFF
from sma_control.safety import SafetyCfg, SafetyMonitor
from sma_control.spring import BiasSpring

TRIP_CODE = {None: 0, "sensor_range": 1, "over_temperature": 2, "tc_disagree": 3, "tc_reads_high": 4, "tc_reads_low": 5, "stall": 6,
             "current_mismatch": 7, "resistance": 8, "watchdog": 9}


def g(x):
    return f"{float(x):.17g}"


def header(sigma_ctl):
    p, cfg = ev.NOM, ev.make_cfg()
    u, M = UniaxialMEFF(p.meff), p.meff
    K = BiasSpring(sigma_ctl, getattr(ev, "K_CTL_OVERRIDE", None) or ev.K_SPRING).K(p)
    sc = SafetyCfg()
    cd = asdict(cfg)
    order = ["Ts", "Kp", "Ki", "k_soft", "dTdt_max", "margin_up", "margin_dn", "T_safe_above_fs", "Kr", "Ki_r", "Kr_dn", "Ki_r_dn", "tau_filt", "Kp_hold",
             "Ki_hold", "drdt_up", "dn_frac", "Kaw", "db_start", "db_done", "T_amb_est", "T_limit_above_Af", "T_abs_max", "T_limit_expand_rate",
             "learn_margin", "hold_step", "hold_drift", "hold_below_As", "sigma_unc", "overshoot_exit",
             "heat_only", "creep_confirm", "creep_rate", "creep_detect", "q_tau", "approach_gain", "approach_span", "approach_min_rate", "fast_detect",
             "creep_pcap", "settle_dr", "settle_s", "full_min", "near_tol", "drop_tol", "reset_r", "reset_still_s", "reset_timeout"]
    so = ["T_trip", "T_disagree", "T_disagree_dyn", "h_min", "h_max", "dT_amb", "band", "T_lat_lo", "stall_dT", "t_stall", "stall_frac", "I_tol_abs",
          "I_tol_rel", "t_I_mismatch", "R_tol", "mc_min_frac", "mc_max_frac", "lat_max_frac", "tc_lag_max"]
    lines = ["wire " + " ".join(g(x) for x in (p.L, p.A, p.surface, p.heat_capacity, p.latent_total, p.rho_eM, p.rho_eA, p.h, p.Imax)),
             "meff " + " ".join(g(x) for x in (M.T0, M.R, M.beta, M.eps_L, M.h, M.m, u.E_A, u.E_M)),
             "cfg " + " ".join(g(float(cd[k])) for k in order),
             "ctl " + g(sigma_ctl) + " " + g(K),
             "safety " + " ".join(g(getattr(sc, k)) for k in so) + f" {g(sc.T_range[0])} {g(sc.T_range[1])} {sc.n_freeze}"]
    return lines


CTL_SEEN = []
_mk = ev.make_ctl


def _capture(*a, **k):
    c = _mk(*a, **k)
    CTL_SEEN.append(c)
    return c


ev.make_ctl = _capture


def export(name, out_dir, fault=None, t_fault=0.0, use_T2=True, sigma_ctl=ev.SIG, **seq_kw):
    cfg = ev.make_cfg()
    guard = ev.Guard(ev.NOM, cfg.Ts, fault=fault, t_fault=t_fault, label=name, use_T2=use_T2,
                     T_lag=getattr(seq_kw.get("sens"), "T_lag", 0.0))
    L, seg, home = ev.sequence(guard=guard, sigma_ctl=sigma_ctl, t_max=seq_kw.pop("t_max", 900.0), **seq_kw)
    T1a, T2a = guard.arm_info
    lines = header(sigma_ctl)
    lines.append(f"arm {g(T1a)} {0 if T2a is None else 1} {g(0.0 if T2a is None else T2a)}")
    lines.append(f"rows {len(guard.rows)}")
    for (phase, tgt, T1, d, I_ctl, I_in, I_meas, T2, V, beat, Ip, relay) in guard.rows:
        lines.append(" ".join([str(phase), g(tgt), g(T1), g(d), g(I_ctl), g(I_in), g(I_meas), "0" if T2 is None else "1", g(0.0 if T2 is None else T2),
                               "1", g(V), "1", str(beat), g(Ip), "1" if relay else "0"]))
    path = Path(out_dir) / f"{name}.txt"
    path.write_text("\n".join(lines) + "\n")
    cal = home.get("cal")
    return path, dict(rows=len(guard.rows), cal=cal, ctl=CTL_SEEN[-1] if CTL_SEEN else None, home_fault=home["fault"], trip=TRIP_CODE.get(guard.mon.reason if guard.mon.tripped else None, -1),
                      reached=len(seg))


def scenarios():
    m = ev.NOM.meff
    return {
        "nominal": dict(),
        "slow_cooling_h5": dict(p_true=replace(ev.NOM, h=5.0), t_max=2500.0),
        "wrong_stress_300": dict(sigma_ctl=300.0),
        "noisy_sensors": dict(sens=ev.Sensors(Tn=0.5, dn=0.1, quant=0.05, T_lag=1.0)),
        "fault_tc_offset": dict(fault="tc_offset", t_fault=290.0),
        "fault_stuck_on": dict(fault="stuck_on", t_fault=290.0),
        "fault_position_frozen": dict(fault="d_stuck", t_fault=250.0),
        "fault_controller_hang": dict(fault="hang", t_fault=250.0),
    }


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "replay")
    out.mkdir(parents=True, exist_ok=True)
    for name, kw in scenarios().items():
        path, info = export(name, out, **kw)
        print(f"{name:24s} {info['rows']:6d} rows  home fault {info['home_fault']}  trip code {info['trip']}  targets {info['reached']}")
