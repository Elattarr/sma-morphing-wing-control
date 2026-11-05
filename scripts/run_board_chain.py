"""Whole chain: Python wire plant  <->  C board layer (slow loop, 20 kHz current loop, wire-voltage estimate, safety monitor) in the PC power-stage simulation.

python scripts/run_board_chain.py nominal|window18|broken_wire|all     prints one summary line per scenario
Run in its own process: the plant is chosen by environment variables read at import (like evaluate_heatonly.py).  Nothing here runs on a board.
"""
import ctypes as C
import os
import subprocess
import sys
import tempfile
from pathlib import Path
os.environ.update(SMA_D_MM="2", SMA_IMAX="35", SMA_SPRING_K="6.36", SMA_SIG="2.2", SMA_NO_HOMING="1")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import evaluate_heatonly as eh
import evaluate_controller as ev

ROOT = Path(__file__).resolve().parents[1]
SLOW = 0.05


def build():
    out = Path(tempfile.mkdtemp()) / "libboard_sim.so"
    c = ROOT / "c"
    subprocess.run(["gcc", "-std=c99", "-O2", "-DSMA_FLOAT", "-shared", "-fPIC", str(c / "board/board.c"), str(c / "board/fast_loop.c"),
                    str(c / "board/hal_sim.c"), str(c / "sma.c"), "-o", str(out), "-lm"], check=True)
    L = C.CDLL(str(out))
    D = C.c_double
    L.sim_set_plant.argtypes = [D, D, D]; L.sim_set_errors.argtypes = [D, D, D, D]; L.sim_set_sensors.argtypes = [D, D, D]
    L.sim_set_faults.argtypes = [C.c_int] * 3; L.sim_slow.argtypes = [D]; L.sim_get.argtypes = [C.POINTER(D)]
    return L


def get(L):
    o = (C.c_double * 20)(); L.sim_get(o); return list(o)


def run_chain(L, t_end=90.0, p_true=None, R_fault_t=None, seed=3):
    p = p_true or ev.NOM
    w = ev.make_plant(p, 22.0, 1.0, ev.SIG)
    rng = np.random.default_rng(seed)
    L.sim_init(); L.sim_set_errors(0.0, 0.0, 0.0, 0.5)
    log, arcs = [], 0
    for k in range(int(t_end / SLOW)):
        t = k * SLOW
        R_w = p.resistance(1.0 - w.xi) * (50.0 if (R_fault_t is not None and t >= R_fault_t) else 1.0)
        L.sim_set_plant(R_w, 5.0, 0.012)
        L.sim_set_sensors(w.T + 0.1 * rng.standard_normal(), w.T + 0.2 * rng.standard_normal(), w.disp_mm(ev.SIG) + 0.02 * rng.standard_normal())
        L.sim_slow(1.0)
        o = get(L)
        ev.plant_step(w, o[1] if o[9] > 0.5 else 0.0, SLOW, ev.SIG, 22.0)
        log.append((t, w.T, 1.0 - w.xi, o[0], o[3], o[4], o[8], o[7], o[1]))
        arcs = o[10]
    return np.array(log), arcs


def summary(name, L_, arcs):
    trips = L_[L_[:, 6] > 0]
    return dict(name=name, rmax=float(L_[:, 2].max()), Tmax=float(L_[:, 1].max()), E=float(np.sum(L_[:, 8] ** 2 * 0.065) * SLOW),
                trip=int(trips[0, 6]) if len(trips) else 0, t_trip=float(trips[0, 0]) if len(trips) else -1.0, arcs=int(arcs), I_end=float(L_[-1, 3]),
                Rest_err_max=float(np.nanmax(np.abs(L_[L_[:, 3] > 5, 5] - 0.0))) if np.any(L_[:, 3] > 5) else 0.0)


def main(which):
    L = build()
    out = []
    if which in ("nominal", "all"):
        l, a = run_chain(L); out.append(summary("nominal", l, a))
        Lp, seg, home = ev.sequence(targets=(1.0,), hold_s=10, t_max=60)
        out[-1]["python_rmax"], out[-1]["python_E"] = float(Lp[:, 2].max()), float(Lp[-1, 6])
    if which in ("window18", "all"):
        l, a = run_chain(L, p_true=eh.wing_window()); out.append(summary("window18", l, a))
    if which in ("broken_wire", "all"):
        l, a = run_chain(L, t_end=40.0, R_fault_t=6.0); out.append(summary("broken_wire", l, a))
    for s in out:
        print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in s.items()))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
