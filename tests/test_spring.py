import numpy as np
from sma_control import SMAParams
from sma_control.cycle_plant import SMAWireCycle
from sma_control.spring import SMAWireSpring, BiasSpring, SpringStress

P = SMAParams()


def heat_cool(w, step_args, T_hot=105.0, t_cool=2500.0, dt=0.05):
    """heat at 15 A until T_hot, then switch off; returns (T, xi) log"""
    log, hot = [], True
    for k in range(int((t_cool + 300) / dt)):
        if hot and w.T >= T_hot:
            hot = False
        w.step(15.0 if hot else 0.0, dt, *step_args)
        log.append((w.T, w.xi))
    return np.array(log)


def test_zero_stiffness_equals_constant_stress_plant():
    a = SMAWireCycle(P, 22.0, 1.0)
    b = SMAWireSpring(P, BiasSpring(250.0, 0.0), 22.0, 1.0)
    dt, hot, mx = 0.05, True, 0.0
    for _ in range(int(2800 / dt)):
        if hot and a.T >= 105.0:
            hot = False
        I = 15.0 if hot else 0.0
        a.step(I, dt, 250.0, 22.0)
        b.step(I, dt, 22.0)
        mx = max(mx, abs(a.T - b.T), abs(a.xi - b.xi))
    assert a.xi > 0.999 and a.T < 22.1                 # the run really went through recovery and re-extension
    assert mx < 1e-6


def test_stress_law_endpoints_and_monotone():
    st = SMAWireSpring(P, BiasSpring(250.0, 20.0)).st
    assert abs(st.sigma(1.0) - 250.0) < 1e-9
    xs = np.linspace(0, 1, 21)
    s = np.array([st.sigma(x) for x in xs])
    assert np.all(np.diff(s) < 0)                        # contracting raises the spring force
    # contraction of the spring: Delta sigma = K * Delta eps  (a few percent of eps_tr correction from the modulus change)
    K = st.K
    assert abs((s[0] - s[-1]) - K * (st.strain(0.0) - st.strain(1.0)) * -1) < 1e-9


def test_strain_inverse():
    st = SMAWireSpring(P, BiasSpring(250.0, 20.0)).st
    for x in np.linspace(0, 1, 11):
        assert abs(st.xi_from_strain(st.strain(x)) - x) < 1e-6


def test_recovery_and_reextension_follow_the_coupled_limits():
    w = SMAWireSpring(P, BiasSpring(250.0, 20.0), 22.0, 1.0)
    log = heat_cool(w, (22.0,))
    T, xi = log[:, 0], log[:, 1]
    heating = np.r_[True, np.diff(xi) < -1e-9][:len(xi)] & (xi > 0.02) & (xi < 0.98)
    cooling = np.r_[False, np.diff(xi) > 1e-9][:len(xi)] & (xi > 0.02) & (xi < 0.98)
    assert heating.sum() > 50 and cooling.sum() > 50
    st = w.st
    assert max(abs(T[i] - st.T_rec(xi[i])) for i in np.where(heating)[0]) < 1e-3
    assert max(abs(T[i] - st.T_fwd(xi[i])) for i in np.where(cooling)[0]) < 1e-3
    assert xi[-1] > 0.999                                  # fully re-extended by the spring at ambient


def test_window_is_wider_than_free_wire_by_the_spring_stress():
    w = SMAWireSpring(P, BiasSpring(250.0, 20.0))
    st, M = w.st, w.p.meff
    width = st.T_rec(0.0) - st.T_rec(1.0)
    expected = M.h * M.eps_L / M.beta + np.sqrt(2 / 3) * (st.sigma(0) - st.sigma(1)) / M.beta
    assert abs(width - expected) < 1e-9
    assert width > 5.0                                     # about 7 K for k = 20 N/mm, against 1.04 K free


def test_spring_preload_threshold_for_return_at_ambient():
    """forward start must stay above 22 C all the way to xi = 1: needs sigma_ext > about 211 MPa (set by the MEFF constants)"""
    for s_ext, returns in ((205.0, False), (225.0, True)):
        w = SMAWireSpring(P, BiasSpring(s_ext, 20.0), 22.0, 1.0)
        log = heat_cool(w, (22.0,), t_cool=3000.0)
        assert (log[-1, 1] > 0.999) == returns, (s_ext, log[-1])
    a = SMAWireSpring(P, BiasSpring(211.2, 20.0)).st.T_fwd(1.0)
    assert abs(a - 22.0) < 0.1
