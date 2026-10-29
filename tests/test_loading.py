import numpy as np
import pytest
from sma_control import SMAParams, RecoveryMap
from sma_control.loading import UniaxialMEFF


def test_onset_and_finish_stress_at_22C():
    u = UniaxialMEFF()
    assert u.sigma_start(22.0) == pytest.approx(203.97, abs=0.05)      # m = 0.01 (204.5 for m = 0)
    assert u.sigma_finish(22.0) - u.sigma_start(22.0) == pytest.approx(7.1, abs=0.05)


def test_loading_reaches_full_transformation_and_unloading_keeps_it():
    u = UniaxialMEFF()
    for s in np.linspace(0, 272.8, 400):
        u.step(s, 22.0)
    assert u.xi == pytest.approx(1.0)
    for s in np.linspace(272.8, 0, 400):
        u.step(s, 22.0)
    assert u.xi == pytest.approx(1.0)                                   # frozen strain after unloading
    assert u.strain(0.0) == pytest.approx(0.03266, abs=2e-5)


def test_sigma_zero_reproduces_recovery_map():
    """regression: at zero stress and xi0 = 1 the stress model is exactly the W1/W2 map"""
    p = SMAParams()
    rm = RecoveryMap(p)
    u = UniaxialMEFF(meff=p.meff, xi=1.0)
    assert u.A_s(0.0) == pytest.approx(p.meff.As)
    assert u.A_f(0.0) == pytest.approx(p.meff.Af)
    for T in np.linspace(40, 60, 81):
        u.step(0.0, T)
        assert 1.0 - u.xi == pytest.approx(rm.r_of_T(T), abs=1e-9)


def test_partial_xi0_shifts_recovery_start_by_one_minus_xi0_times_window():
    p = SMAParams()
    u = UniaxialMEFF(meff=p.meff, xi=0.954)
    start = None
    for T in np.linspace(50.0, 52.0, 4001):
        u.step(0.0, T)
        if u.xi < 0.954 - 1e-9 and start is None:
            start = T
    assert start - p.meff.As == pytest.approx(0.046 * p.meff.window, abs=1e-3)


def test_stress_shifts_transformation_temperatures():
    u = UniaxialMEFF()
    assert u.shift_per_MPa == pytest.approx(0.1458, abs=1e-3)
    assert u.A_s(100.0) - u.A_s(0.0) == pytest.approx(14.58, abs=0.05)
    assert u.A_f(100.0) - u.A_s(100.0) == pytest.approx(u.A_f(0.0) - u.A_s(0.0))   # window does not change at constant stress
