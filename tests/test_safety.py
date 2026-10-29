"""Independent safety monitor: each check on its own, latch behaviour, and no false trip on a healthy run."""
import numpy as np
import pytest
from sma_control import SMAParams
from sma_control.safety import SafetyMonitor, SafetyCfg

P = SMAParams(d=0.5e-3, Imax=3.0)
TS = 0.05


def mon(**kw):
    m = SafetyMonitor(P, TS, SafetyCfg(**kw))
    m.arm(22.0, 22.0)
    return m


def run(m, n, **kw):
    out = None
    for _ in range(n):
        args = dict(I_cmd=0.0, I_meas=0.0, T1=22.0, T2=22.0, d=0.0, V=None, beat=None)
        args.update({k: (v(m.t) if callable(v) else v) for k, v in kw.items()})
        out = m.step(**args)
    return out


def test_healthy_idle_does_not_trip():
    m = mon()
    assert run(m, 400, beat=lambda t: int(t / TS)) == (0.0, True) and not m.tripped


def test_arm_rejects_open_thermocouple_rail():
    m = SafetyMonitor(P, TS)
    m.arm(-40.0, 22.0)
    assert m.tripped and m.reason == "sensor_range"


def ramp(T_end, rate=30.0):
    """physically consistent heating signal: 3 A on a 0.5 mm wire heats at about 35 K/s"""
    return lambda t: min(T_end, 22.0 + rate * t)


def test_hard_over_temperature_on_either_channel():
    m = mon()
    run(m, 120, I_cmd=3.0, I_meas=3.0, T1=ramp(150.0), T2=ramp(150.0))
    assert m.tripped and m.reason == "over_temperature"
    m = mon()                                          # second channel alone over the limit, first reads plausible
    run(m, 120, I_cmd=3.0, I_meas=3.0, T1=ramp(100.0), T2=ramp(150.0))
    assert m.tripped and m.reason in ("over_temperature", "tc_disagree")


def test_trip_latches_and_reset_clears():
    m = mon()
    run(m, 1, T1=140.0, T2=140.0)
    assert run(m, 5, I_cmd=2.0) == (0.0, False)
    m.reset()
    assert not m.tripped


def test_channel_disagreement():
    m = mon()
    run(m, 5, T1=22.0, T2=40.0)
    assert m.reason == "tc_disagree"


def test_envelope_catches_open_thermocouple_with_power_applied():
    """single sensor reading 0 C while 3 A flows: the lower bound rises above it within a few seconds"""
    m = mon()
    run(m, 400, I_cmd=3.0, I_meas=3.0, T1=0.0, T2=None)
    assert m.tripped and m.reason == "tc_reads_low"


def test_current_mismatch_power_stage_stuck_on():
    m = mon()
    run(m, 100, I_cmd=0.0, I_meas=3.0)
    assert m.reason == "current_mismatch"


def test_wire_resistance_out_of_range():
    m = mon()
    run(m, 5, I_cmd=1.0, I_meas=1.0, V=3.0 * P.resistance(0.5) * 1.0)
    run(m, 1, I_cmd=1.0, I_meas=1.0, V=3.0 * P.resistance(0.5) * 1.0)
    assert m.reason == "resistance"


def test_watchdog_on_frozen_heartbeat():
    m = mon()
    run(m, 5, beat=lambda t: int(t / TS))
    run(m, 20, beat=7)
    assert m.reason == "watchdog"


def test_stall_needs_no_motion_above_af_and_resets_on_motion():
    kw = dict(I_cmd=2.0, I_meas=2.0, T1=ramp(100.0, 8.0), T2=ramp(100.0, 8.0))   # 2 A heats at about 15 K/s
    m = mon()
    m.set_calibration(stroke_mm=6.0, Af_max=94.0)
    run(m, 600, d=lambda t: 0.05 * t, **kw)             # moving 0.05 mm/s = 1 mm per 20 s: more than 15 % of 6 mm per 8 s? no, 0.4 mm
    assert m.reason == "stall"
    m = mon()
    m.set_calibration(stroke_mm=6.0, Af_max=94.0)
    run(m, 600, d=lambda t: 0.3 * t, **kw)              # 2.4 mm per 8 s: real motion keeps resetting the timer
    assert m.reason != "stall"


def test_no_false_trip_in_healthy_homing_and_sequence():
    import os, sys
    os.environ.update(SMA_D_MM="0.5", SMA_IMAX="3", SMA_SPRING_K="1.25")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import evaluate_controller as ev
    L, seg, home = ev.sequence(label="healthy")
    assert len(seg) == 4 and not home["guard"].mon.tripped
