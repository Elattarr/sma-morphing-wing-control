"""Homing / calibration tests on the 0.5 mm spring-biased wire (the evaluation harness builds plant, sensors and controller)."""
import os, sys
os.environ.update(SMA_D_MM="0.5", SMA_IMAX="3", SMA_SPRING_K="1.25")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from dataclasses import replace
import pytest
import evaluate_controller as ev


def _home(p_true=None, sigma_true=ev.SIG, sigma_ctl=ev.SIG, T_amb=22.0):
    cfg = ev.make_cfg()
    w = ev.make_plant(p_true or ev.NOM, T_amb, 1.0, sigma_true)
    ctl = ev.make_ctl(sigma_ctl, cfg)
    home = ev.do_homing(w, ctl, ev.Sensors(T0=T_amb), cfg, T_amb, sigma_true)
    return w, home


def test_homing_measures_windows_of_the_true_wire():
    w, home = _home()
    cal = home["cal"]
    assert home["fault"] is None
    assert cal.As == pytest.approx(w.st.T_rec(1.0), abs=2.0)
    assert cal.Af == pytest.approx(w.st.T_rec(0.0), abs=2.5)
    assert cal.Tfs == pytest.approx(w.st.T_fwd(0.0), abs=2.0)


def test_homing_slow_cooling_waits_for_the_extension_tail():
    """h_c = 5: the extension tail moves less than the threshold per 2 s, homing must still wait for the real end"""
    w, home = _home(replace(ev.NOM, h=5.0))
    assert home["fault"] is None
    assert w.xi > 0.98                                   # wire is back at the extended end (xi = 1)


def test_homing_refuses_when_preload_cannot_return_the_wire():
    m = ev.NOM.meff
    wide = replace(ev.NOM, meff=replace(m, h=12.0 * m.beta / m.eps_L))
    _, home = _home(wide)                               # 250 MPa preload: T_fwd(1) is below ambient
    assert home["fault"] == "preload_too_low"
    assert home["cal"] is None


@pytest.mark.parametrize("sig_ctl", [200.0, 300.0])
def test_wrong_stress_assumption_is_corrected_by_homing(sig_ctl):
    L, seg, home = ev.sequence(sigma_ctl=sig_ctl)
    assert home["fault"] is None and len(seg) == 4
