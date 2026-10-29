import numpy as np
from sma_control import SMAParams
from sma_control.cycle_plant import SMAWireCycle
from sma_control.two_way import SMAWireTwoWay, TwoWayParams


def run_cool(w, sigma, t=900.0, dt=0.05):
    xs = []
    for _ in range(int(t / dt)):
        w.step(0.0, dt, sigma, 22.0)
        xs.append((w.T, w.xi))
    return np.array(xs)


def heat(w, sigma, I=15.0, dt=0.02):
    Toff = w.u.A_f(sigma) + 3.0
    while w.T < Toff:
        w.step(I, dt, sigma, 22.0)


def test_two_way_cooling_returns_capped_fraction_at_zero_stress():
    w = SMAWireTwoWay(SMAParams(), T0=22.0, xi0=1.0, tw=TwoWayParams(Ms=35.0, xi_cap=0.3))
    heat(w, 0.0)
    xs = run_cool(w, 0.0)
    i = int(np.argmax(xs[:, 1] > 1e-3))
    assert abs(xs[i, 0] - 35.0) < 0.15
    assert abs(xs[-1, 1] - 0.3) < 1e-6


def test_two_way_does_not_change_recovery_branch():
    a = SMAWireCycle(SMAParams(), T0=35.0, xi0=1.0)
    b = SMAWireTwoWay(SMAParams(), T0=35.0, xi0=1.0)
    for _ in range(300):
        a.step(15.0, 0.01, 0.0, 22.0); b.step(15.0, 0.01, 0.0, 22.0)
    assert abs(a.T - b.T) < 1e-9 and abs(a.xi - b.xi) < 1e-9
