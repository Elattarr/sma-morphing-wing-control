"""Heat-only (one-way) controller: behaviour tests on the spring plant (2 mm thesis wire, structure stiffness only)."""
import numpy as np
import pytest
from dataclasses import replace
from sma_control.params import SMAParams
from sma_control.spring import SMAWireSpring, BiasSpring
from sma_control.elements import SMAWireElements
from sma_control.cycle_controller import CycleController, cfg_heat_only
from sma_control.safety import SafetyMonitor
from sma_control.te_map import TEMap
from sma_control.hardware_io import RigIO, LoopState, run_step

P = SMAParams(d=2e-3, Imax=35.0)
SIG, KS, TS = 2.2, 6.36, 0.05


class Rig:
    def __init__(self, p_true=P, struct=None, T0=22.0, sensor_d=None):
        self.ks = BiasSpring(SIG, KS)
        self.w = SMAWireElements(p_true, self.ks, T0=T0, xi0=1.0, **struct) if struct else SMAWireSpring(p_true, self.ks, T0=T0, xi0=1.0)
        self.ctl = CycleController(P, SIG, cfg_heat_only(P, Ts=TS), K_spring=self.ks.K(P))
        self.t, self.log = 0.0, []

    def step(self, tgt, n=1, T_amb=22.0):
        for _ in range(n):
            I = self.ctl.step(tgt, self.w.T, self.w.disp_mm(SIG))
            self.w.step(I, TS, T_amb, 1.0)
            self.t += TS
            self.log.append((self.t, self.w.T, 1.0 - self.w.xi, tgt, I))
        return self.log[-1]

    def r(self):
        return 1.0 - self.w.xi


def test_moves_up_to_a_target_and_holds_with_zero_current():
    g = Rig()
    g.step(0.5, int(40 / TS))
    assert abs(g.r() - 0.5) < 0.02 and g.ctl.mode == "hold"
    assert g.log[-1][4] == 0.0                                          # hold = power off
    assert max(x[1] for x in g.log) < 60.0


def test_never_enters_down_mode_when_target_is_lowered():
    g = Rig()
    g.step(0.6, int(40 / TS))
    r_before = g.r()
    modes = set()
    for _ in range(int(30 / TS)):
        g.step(0.2)
        modes.add(g.ctl.mode)
    assert "down" not in modes and abs(g.r() - r_before) < 0.01
    assert g.ctl.fault is None and max(x[4] for x in g.log[-int(30 / TS):]) == 0.0


def test_reaches_a_higher_target_after_a_lower_one_without_overshoot():
    g = Rig()
    for tgt in (0.25, 0.5, 0.75, 1.0):
        g.step(tgt, int(25 / TS))
        assert abs(g.r() - tgt) < 0.02
        assert max(x[2] for x in g.log) - tgt < 0.02


def test_accepted_full_stroke_is_not_chased_again():
    """r settles just under 1.0: the controller must not creep towards the temperature cap trying to reach the last 1 %"""
    g = Rig()
    g.step(1.0, int(120 / TS))
    assert g.ctl.mode == "hold" and max(x[1] for x in g.log) < 66.0 and g.ctl.fault is None


def test_wide_window_is_followed_and_stopped_at_the_end_of_travel():
    g = Rig(struct=dict(spread_K=4.0))
    g.step(1.0, int(120 / TS))
    assert g.r() > 0.97 and max(x[1] for x in g.log) < 66.0 and g.ctl.mode == "hold"


def test_first_stroke_learner_on_a_short_stroke_wire():
    M = P.meff
    g = Rig(p_true=replace(P, meff=replace(M, eps_L=M.eps_L * 0.9)))
    g.step(0.5, int(30 / TS))
    g.step(1.0, int(80 / TS))
    c = g.ctl
    assert c.cal is not None and c.saturated and c.fault is None
    assert 52.0 < c.cal.Af < 66.0 and c.cal.eps_ext != c.cal.eps_con
    assert c.r_est > 0.95                                                # the estimate is on the learned mapping


def test_far_too_short_stroke_is_a_fault_not_a_hang():
    M = P.meff
    g = Rig(p_true=replace(P, meff=replace(M, eps_L=M.eps_L * 0.6)))
    g.step(1.0, int(100 / TS))
    assert g.ctl.fault == "stalled_short" and g.log[-1][4] == 0.0


def test_unexpected_drop_is_a_fault_and_power_stays_off():
    g = Rig()
    g.step(0.6, int(40 / TS))
    g.w.xi = 1.0                                                        # somebody pushed the wing back
    g.step(0.6, int(5 / TS))
    assert g.ctl.fault == "position_dropped"
    assert all(x[4] == 0.0 for x in g.log[-int(2 / TS):])


def test_reset_rearms_after_the_wing_is_back_and_times_out_otherwise():
    g = Rig()
    g.step(0.6, int(40 / TS))
    g.w.xi = 1.0
    g.step(0.6, int(3 / TS))
    assert g.ctl.fault == "position_dropped"
    g.ctl.begin_reset()
    assert g.ctl.resetting and not g.ctl.ready
    g.step(0.6, int(10 / TS))                                           # wire is back at r = 0: re-armed after reset_still_s
    assert g.ctl.ready and g.ctl.fault is None
    g.step(0.6, int(40 / TS))
    assert abs(g.r() - 0.6) < 0.02
    # timeout: wing never returns
    g2 = Rig()
    g2.ctl.c.reset_timeout = 5.0
    g2.step(0.6, int(40 / TS))
    g2.w.xi = 0.4                                                       # held at r = 0.6
    g2.ctl.fault = "position_dropped"
    g2.ctl.begin_reset()
    g2.step(0.6, int(8 / TS))
    assert g2.ctl.fault == "reset_timeout"


def test_motion_detection_ignores_noise_spikes_before_the_window():
    """0.02 mm displacement noise (0.3 % of the stroke), cold wire: no false motion in 10 seeds.
    Measured limit of the current detector: 0.03 mm -> 1/20 false detections, 0.05 mm -> 4/20 (see STATUS)."""
    for seed in range(10):
        g = Rig()
        rng = np.random.default_rng(seed)
        for _ in range(int(3 / TS)):
            d = g.w.disp_mm(SIG) + 0.02 * rng.standard_normal()
            I = g.ctl.step(0.5, g.w.T, d)
            g.w.step(I, TS, 22.0, 1.0)
        assert g.ctl.T_start_obs is None, seed


def test_te_map_roundtrip_and_validation():
    m = TEMap()
    assert m.r_to_te(0.5) == pytest.approx(-21.1) and m.te_to_r(-31.65) == pytest.approx(0.75)
    m2 = TEMap(((0.0, 0.0), (0.5, -30.0), (1.0, -42.2)))
    for r in (0.0, 0.2, 0.5, 0.9, 1.0):
        assert m2.te_to_r(m2.r_to_te(r)) == pytest.approx(r)
    with pytest.raises(ValueError):
        TEMap(((0.0, 0.0), (0.5, -10.0), (0.5, -20.0), (1.0, -42.2)))


class SimRig(RigIO):
    def __init__(self, g):
        self.g, self.I, self.relay = g, 0.0, True
    def read_T1(self): return self.g.w.T
    def read_T2(self): return self.g.w.T
    def read_disp_mm(self): return self.g.w.disp_mm(SIG)
    def read_I_meas(self): return self.I if self.relay else 0.0
    def write_current(self, amps): self.I = amps
    def write_relay(self, closed): self.relay = closed


def test_io_loop_runs_a_move_in_millimetres_with_the_monitor():
    g = Rig()
    io, mon, st, te = SimRig(g), SafetyMonitor(P, TS), LoopState(), TEMap()
    mon.arm(22.0, 22.0)
    for _ in range(int(40 / TS)):
        run_step(io, g.ctl, mon, -21.1, te, st)
        g.w.step(io.I, TS, 22.0, 1.0)
    assert abs(te.r_to_te(g.r()) - (-21.1)) < 1.0 and not mon.tripped and st.fault is None


def test_safety_catches_a_frozen_thermocouple_and_not_a_noisy_one():
    from sma_control.safety import SafetyCfg
    mon = SafetyMonitor(P, TS)
    mon.arm(22.0, None)
    rng = np.random.default_rng(0)
    for k in range(int(30 / TS)):                         # heated, noisy reading: no trip
        mon.step(5.0, 5.0, 40.0 + 0.1 * rng.standard_normal())
    assert not mon.tripped or mon.reason != "tc_frozen"
    mon2 = SafetyMonitor(P, TS)
    mon2.arm(22.0, None)
    for k in range(int(30 / TS)):                         # the same value for ever while heated
        mon2.step(5.0, 5.0, 30.0)
        if mon2.tripped:
            break
    assert mon2.tripped and mon2.reason in ("tc_frozen", "tc_reads_low")
    mon3 = SafetyMonitor(P, TS)
    mon3.arm(22.0, None)
    for k in range(int(60 / TS)):                         # cold idle wire at ambient, heater off: constant reading is fine
        mon3.step(0.0, 0.0, 22.0)
    assert not mon3.tripped
