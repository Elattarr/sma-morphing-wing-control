"""Replay scenarios of the heat-only controller on the 2 mm thesis wire (no homing) for the C port; run in its own process
(the plant is chosen by environment variables read at import).

python scripts/export_replay_heatonly.py <out_dir>
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_heatonly as h            # sets the 2 mm plant and patches the harness
import export_replay as ex
from dataclasses import replace
import evaluate_controller as ev

UP = h.UP
SCEN = {
    "ho_nominal": dict(targets=UP, hold_s=10.0, t_max=120.0),
    "ho_spread4": dict(targets=UP, hold_s=10.0, t_max=120.0, struct=dict(spread_K=4.0)),
    "ho_noisy": dict(targets=UP, hold_s=10.0, t_max=120.0, sens=ev.Sensors(Tn=0.3, dn=0.01, T_lag=0.5)),
    "ho_short_stroke": dict(targets=(1.0,), hold_s=10.0, t_max=150.0,
                            p_true=replace(ev.NOM, meff=replace(ev.NOM.meff, eps_L=ev.NOM.meff.eps_L * 0.9))),
    "ho_fault_tc_offset": dict(targets=UP, hold_s=10.0, t_max=120.0, fault="tc_offset", t_fault=20.0),
    "ho_fault_stuck_on": dict(targets=UP, hold_s=10.0, t_max=120.0, fault="stuck_on", t_fault=20.0),
}

if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "replay_ho")
    out.mkdir(parents=True, exist_ok=True)
    for name, kw in SCEN.items():
        path, info = ex.export(name, out, **kw)
        c = info["ctl"]
        print(f"{name:22s} {info['rows']:6d} rows  trip {info['trip']}  targets {info['reached']}  py_cal {c.cal is not None} py_fault {c.fault}")
