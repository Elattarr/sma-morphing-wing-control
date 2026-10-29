import numpy as np
from sma_control import SMAParams, SMAWire
from sma_control.cycle_plant import SMAWireCycle


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


def test_heating_matches_one_way_plant_at_zero_stress():
    p = SMAParams()
    a, b = SMAWireCycle(p, T0=35.0, xi0=1.0), SMAWire(p, T0=35.0)
    for _ in range(int(6.0 / 0.01)):
        a.step(15.0, 0.01, 0.0, 22.0)
        b.step(15.0, 0.01, 22.0)
    assert abs(a.T - b.T) < 0.05
    assert abs((1.0 - a.xi) - b.r) < 0.01


def test_one_way_cooling_at_zero_stress_changes_nothing():
    w = SMAWireCycle(SMAParams(), T0=22.0, xi0=1.0)
    heat(w, 0.0)
    assert w.xi < 1e-6
    xs = run_cool(w, 0.0)
    assert xs[-1, 0] < 22.1 and xs[:, 1].max() < 1e-6


def test_bias_stress_forward_start_matches_theory():
    for sigma in (250.0, 300.0):
        w = SMAWireCycle(SMAParams(), T0=22.0, xi0=1.0)
        heat(w, sigma)
        xs = run_cool(w, sigma)
        i = int(np.argmax(xs[:, 1] > 1e-3))
        assert xs[i, 1] > 1e-3
        assert abs(xs[i, 0] - w.T_forward_start(sigma)) < 0.15
        assert xs[-1, 1] > 0.999                              # fully re-extended at 22 C


def test_bias_below_threshold_never_returns_at_ambient():
    w = SMAWireCycle(SMAParams(), T0=22.0, xi0=1.0)
    heat(w, 100.0)
    xs = run_cool(w, 100.0)
    assert xs[:, 1].max() < 1e-6                             # forward start at 6.8 C is below ambient
