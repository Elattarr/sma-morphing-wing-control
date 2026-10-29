import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from sma_control import SMAParams
from sma_control.cycle_plant import SMAWireCycle
from sma_control.cycle_controller import CycleController
import control_bias_cycle as cbc


def test_position_estimate_inverts_the_plant():
    p = SMAParams()
    ctl = CycleController(p, 250.0)
    w = SMAWireCycle(p, T0=22.0, xi0=1.0)
    for xi in (0.0, 0.3, 0.7, 1.0):
        w.xi = xi
        assert abs(ctl.r_from_disp(w.disp_mm(250.0)) - (1.0 - xi)) < 1e-6


def test_window_temperatures_match_the_plant():
    ctl = CycleController(SMAParams(), 250.0)
    w = SMAWireCycle(SMAParams())
    assert abs(ctl.Tfs - w.T_forward_start(250.0)) < 1e-9
    assert abs(ctl.As - w.u.A_s(250.0)) < 1e-9


def test_contract_then_extend_to_half_under_250MPa():
    L, seg, ctl = cbc.run(250.0, 250.0, targets=[1.0, 0.5], t_max=1500.0)
    assert len(seg) == 2
    for tgt, tr, te, E in seg:
        m = (L[:, 0] >= tr) & (L[:, 0] <= te)
        assert abs(L[m, 2].mean() - tgt) < 0.02
    assert L[:, 1].max() < 90.0                                      # never above A_f + 1.6 K


# ---------------------------------------------------------------- bias spring
from sma_control.spring import SMAWireSpring, BiasSpring


def test_spring_position_estimate_inverts_the_plant():
    p = SMAParams()
    sp = BiasSpring(250.0, 20.0)
    ctl = CycleController(p, 250.0, K_spring=sp.K(p))
    w = SMAWireSpring(p, sp)
    for xi in (0.0, 0.3, 0.7, 1.0):
        w.xi = xi
        assert abs(ctl.r_from_disp(w.disp_mm()) - (1.0 - xi)) < 1e-6


def test_spring_window_temperatures_match_the_plant_limits():
    p = SMAParams()
    sp = BiasSpring(250.0, 20.0)
    ctl = CycleController(p, 250.0, K_spring=sp.K(p))
    w = SMAWireSpring(p, sp)
    for xi in (0.0, 1.0):
        w.xi = xi
        s = w.sigma
        assert abs(w.u.xi_r(s, (ctl.As, ctl.Af)[1 - int(xi)]) - xi) < 1e-9     # recovery passes xi at T_rec(xi)
    assert abs(w.u.xi_f(w.st.sigma(0.0), ctl.Tfs)) < 1e-9
    assert ctl.Af - ctl.As > 5.0                                                 # wider than the free-wire 1 K


def test_spring_contract_then_extend_to_half():
    L, seg, ctl = cbc.run(250.0, 250.0, targets=[1.0, 0.5], t_max=1500.0, k_true=20.0)
    assert len(seg) == 2
    for tgt, tr, te, E in seg:
        m = (L[:, 0] >= tr) & (L[:, 0] <= te)
        assert abs(L[m, 2].mean() - tgt) < 0.02
    assert L[:, 1].max() < ctl.Af + 3.0


# ---------------------------------------------------------------- gain scaling with wire size
from dataclasses import replace
from sma_control.cycle_controller import cfg_for_wire, CycleCfg


def test_gains_scale_with_heat_capacity_and_cooling_time():
    ref = SMAParams()
    thin = replace(ref, d=0.5e-3)
    a, b, base = cfg_for_wire(ref), cfg_for_wire(thin), CycleCfg()
    assert abs(a.Kp - base.Kp) < 1e-12 and abs(a.tau_filt - base.tau_filt) < 1e-12        # reference wire unchanged
    k = thin.heat_capacity / ref.heat_capacity
    assert abs(b.Kp - base.Kp * k) < 1e-12 and abs(b.Ki_hold - base.Ki_hold * k) < 1e-12
    assert abs(thin.heat_capacity / b.Kp - ref.heat_capacity / base.Kp) < 1e-9            # same loop time constant m c / Kp
    tr = thin.cooling_tau / ref.cooling_tau
    assert abs(b.tau_filt - base.tau_filt * tr) < 1e-12 and abs(b.Kr_dn - base.Kr_dn / tr) < 1e-12
