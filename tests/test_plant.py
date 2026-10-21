from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from sma_control import SMAParams, SMAWire

ROOT = Path(__file__).resolve().parents[1]


def heat_to(T_target, I, p=None, dt=0.001, t_max=200):
    w = SMAWire(p or SMAParams(Imax=50, noise_mm=0))
    t = 0.0
    while w.T < T_target and t < t_max:
        w.step(I, dt); t += dt
    return w, t


def test_window_from_meff_constants():
    m = SMAParams().meff
    assert m.As == pytest.approx(50.79, abs=0.01)
    assert m.Af == pytest.approx(51.82, abs=0.01)
    assert m.window == pytest.approx(1.04, abs=0.01)


def test_map_zero_below_as_one_above_af_and_monotonic():
    w = SMAWire(SMAParams())
    T = np.linspace(20, 70, 500)
    r = np.array([w.map.r_of_T(t) for t in T])
    assert r[0] == 0.0 and r[-1] == 1.0
    assert np.all(np.diff(r) >= 0)


def test_frozen_below_as_and_full_stroke_above_af():
    w, _ = heat_to(45.0, 20)
    assert w.stroke_mm() == 0.0
    w, _ = heat_to(60.0, 20)
    assert w.stroke_mm() == pytest.approx(6.64, abs=0.01)


def test_one_way_no_return_on_cooling():
    w, _ = heat_to(51.3, 20)                  # inside the window
    y = w.stroke_mm()
    assert 0.0 < y < 6.64
    for _ in range(60000):
        w.step(0.0, 0.001)                    # cool for 60 s
    assert w.T < 30
    assert w.stroke_mm() == pytest.approx(y)


def test_energy_balance_closes():
    p = SMAParams(Imax=50, noise_mm=0)
    w = SMAWire(p)
    dt, loss = 0.001, 0.0
    for _ in range(4000):
        loss += p.h * p.surface * (w.T - p.T_inf) * dt
        w.step(25.0, dt)
    stored = p.heat_capacity * (w.T - p.T_inf) + p.latent_total * w.r
    assert w.E_joule == pytest.approx(stored + loss, rel=0.01)


def test_map_matches_ansys_within_0p25K():
    d = pd.read_csv(ROOT / "data" / "ansys_3d_wire_merged.csv")
    s = d[d.step.isin([4, 5]) & (d["T"] > 49)].sort_values("T")
    r = 1.0 - s.p_max.values / s.p_max.values.max()
    inside = (r > 0.05) & (r < 0.95)
    w = SMAWire(SMAParams())
    dT = np.array([w.map.T_of_r(x) for x in r[inside]]) - s["T"].values[inside]
    assert np.abs(dT).max() < 0.25


@pytest.mark.parametrize("I, E_thesis, tol", [(35, 173, 0.03), (26, 182, 0.04), (19, 199, 0.06)])
def test_heating_energy_close_to_thesis_cases(I, E_thesis, tol):
    w, _ = heat_to(73.0, I)
    assert w.E_joule == pytest.approx(E_thesis, rel=tol)
