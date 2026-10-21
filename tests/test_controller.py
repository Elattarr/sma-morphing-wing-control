import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from step_test import run
from sma_control import SMAParams, SMAWire
from sma_control.controller import Controller


@pytest.mark.parametrize("rt", [0.3, 0.5, 1.0])
def test_nominal_step_no_overshoot_and_accurate(rt):
    log, w = run(rt)
    r = log[:, 2]
    assert r.max() - rt <= 0.02                      # <= 2 % of full stroke (0.13 mm)
    assert abs(r[-1] - rt) <= 0.015                  # R1: 0.1 mm = 1.5 % of 6.64 mm
    assert log[:, 1].max() < 56.0                    # R4


@pytest.mark.parametrize("kw", [dict(h_scale=0.2), dict(T_amb=35.0), dict(T_amb=15.0)])
def test_half_step_robust_to_disturbances(kw):
    log, w = run(0.5, **kw)
    assert log[:, 2].max() - 0.5 <= 0.03
    assert abs(log[-1, 2] - 0.5) <= 0.03


@pytest.mark.parametrize("shift", [-1.0, 1.0])
def test_robust_to_map_shift(shift):
    w = SMAWire(SMAParams(T_shift=shift), seed=1)
    ctl = Controller(SMAParams())
    rng = np.random.default_rng(5)
    for _ in range(3000):
        w.step(ctl.step(0.5, w.T + 0.1 * rng.standard_normal(), w.measure()), ctl.c.Ts)
    assert abs(w.r - 0.5) <= 0.03
