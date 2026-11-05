"""Board layer (c/board): 20 kHz current loop, wire-voltage estimate, relay sequencing and the whole chain against the Python wire.
The power stage and the sensors are the PC simulation in c/board/hal_sim.c.  Nothing here runs on a real board."""
import ctypes as C
import os
import shutil
import subprocess
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which("gcc") is None, reason="needs gcc")
SLOW = 0.05
SIM_TRIP_RESISTANCE = 8


@pytest.fixture(scope="module")
def lib(tmp_path_factory):
    out = tmp_path_factory.mktemp("board") / "libboard_sim.so"
    c = ROOT / "c"
    subprocess.run(["gcc", "-std=c99", "-Wall", "-Wextra", "-pedantic", "-Werror", "-O2", "-DSMA_FLOAT", "-shared", "-fPIC",
                    str(c / "board/board.c"), str(c / "board/fast_loop.c"), str(c / "board/hal_sim.c"), str(c / "sma.c"), "-o", str(out), "-lm"], check=True)
    L = C.CDLL(str(out))
    D = C.c_double
    L.sim_set_plant.argtypes = [D, D, D]
    L.sim_set_errors.argtypes = [D, D, D, D]
    L.sim_set_sensors.argtypes = [D, D, D]
    L.sim_set_faults.argtypes = [C.c_int] * 3
    L.sim_slow.argtypes = [D]
    L.sim_direct.argtypes = [D, D, C.c_int, C.c_int, C.POINTER(D)]
    L.sim_get.argtypes = [C.POINTER(D)]
    L.sim_uart.restype = C.c_char_p
    return L


def get(lib):
    o = (C.c_double * 20)()
    lib.sim_get(o)
    return list(o)


def direct(lib, I, R_wire, R_ff, Vbus=5.0, n=400):
    lib.sim_init()
    lib.sim_set_plant(R_wire, Vbus, 0.012)
    tr = (C.c_double * n)()
    lib.sim_direct(I, R_ff, 1, n, tr)
    return np.array(tr)


@pytest.mark.parametrize("R_wire", [0.05, 0.065, 0.09])
@pytest.mark.parametrize("R_ff", [0.062, 0.08, 0.12])
@pytest.mark.parametrize("Vbus", [4.5, 5.0, 5.5])
def test_current_step(lib, R_wire, R_ff, Vbus):
    """35 A step with a wrong feed-forward resistance (0.062 to 0.12 ohm for a wire between 0.062 and 0.102 ohm total) and a bus between 4.5 and 5.5 V"""
    tr = direct(lib, 35.0, R_wire, R_ff, Vbus)
    if Vbus * 0.95 / (R_wire + 0.012) < 35.0:
        pytest.skip("bus too low for 35 A at the maximum duty")
    assert tr.max() < 39.0                                  # below the 40 A software latch
    i = np.where(abs(tr - 35.0) < 0.35)[0]
    assert len(i) and i[0] * 50e-6 < 4e-3                   # settled within 1 % in 4 ms (the slow loop is 50 ms)
    assert abs(tr[-50:].mean() - 35.0) < 0.2
    assert get(lib)[11] == 0


def test_resolution_is_16_mA(lib):
    """one ADC count is 3.3 / 4095 / (50 x 1 mohm) = 16.1 mA"""
    tr = direct(lib, 10.0, 0.06, 0.072)
    assert abs(tr[-100:].mean() - 10.0) < 0.05


def test_overcurrent_latch(lib):
    """a shorted wire: the loop would drive the full duty into 0.022 ohm (216 A); the software latch must cut it within a few samples"""
    tr = direct(lib, 35.0, 0.010, 0.5, n=200)
    assert get(lib)[11] == 1
    assert tr.max() < 80.0 and tr[-1] < 1.0


def test_voltage_estimate_without_a_voltage_sensor(lib):
    """V_wire = mean(D Vbus) - I R_other has to give the wire resistance to a few percent without any extra sensor"""
    for R_wire in (0.055, 0.065, 0.08):
        lib.sim_init(); lib.sim_set_plant(R_wire, 5.0, 0.012); lib.sim_set_errors(0.0, 0.0, 0.0, 0.5)
        lib.sim_set_sensors(22.0, 22.0, 6.5415)
        for _ in range(3):
            lib.sim_slow(0.0)                                # arm
        # drive the current loop directly for a slow period, then let the board take its sample
        lib.sim_direct(20.0, R_wire + 0.012, 1, 3000, None)
        lib.sim_slow(0.0)
        o = get(lib)
        assert abs(o[4] - R_wire) / R_wire < 0.10, (R_wire, o[4])


def test_voltage_estimate_error_budget(lib):
    """bad R_other (+/-50 %) and 1 % bus-divider / 1 % INA gain error: the resistance error at 35 A stays below the 45 % check tolerance"""
    worst = 0.0
    for R_other_true in (0.006, 0.018):
        for g in (-0.01, 0.01):
            lib.sim_init(); lib.sim_set_plant(0.065, 5.0, R_other_true); lib.sim_set_errors(g, 0.0, -g, 0.5)
            lib.sim_set_sensors(22.0, 22.0, 6.5415)
            lib.sim_slow(0.0); lib.sim_slow(0.0)
            lib.sim_direct(33.0, 0.077, 1, 3000, None)
            lib.sim_slow(0.0)
            worst = max(worst, abs(get(lib)[4] - 0.065) / 0.065)
    assert worst < 0.25, worst


# ------------------------------------------------------------------ whole chain: Python wire + C board (own process, plant chosen by environment)
@pytest.fixture(scope="module")
def chain():
    out = subprocess.run([sys.executable, str(ROOT / "scripts/run_board_chain.py"), "all"], capture_output=True, text=True, timeout=600,
                         env=dict(os.environ, PYTHONPATH=str(ROOT / "src")))
    assert out.returncode == 0, out.stderr[-2000:]
    res = {}
    for line in out.stdout.strip().splitlines():
        kv = dict(x.split("=") for x in line.split())
        res[kv["name"]] = {k: float(v) for k, v in kv.items() if k != "name"}
    return res


def test_chain_reaches_full_stroke_like_the_python_controller(chain):
    r = chain["nominal"]
    assert r["arcs"] == 0 and r["trip"] == 0
    assert r["rmax"] > 0.95 and r["Tmax"] < 60.0
    assert abs(r["rmax"] - r["python_rmax"]) < 0.04
    assert 0.6 < r["E"] / r["python_E"] < 1.4


def test_chain_wire_window_18K(chain):
    r = chain["window18"]
    assert r["arcs"] == 0 and r["trip"] == 0                # the stall check follows the observed motion (stall_learn_frac)
    assert r["rmax"] > 0.97 and r["Tmax"] < 73.0


def test_chain_broken_wire_trips_and_opens_relay_without_arc(chain):
    r = chain["broken_wire"]
    assert r["trip"] in (7, 8)                              # resistance (8) or current_mismatch (7)
    assert r["t_trip"] - 6.0 < 2.0
    assert r["arcs"] == 0 and r["I_end"] < 0.5


def test_chain_thermocouple_open_trips(lib):
    lib.sim_init()
    lib.sim_set_sensors(22.0, 22.0, 6.5415)
    for _ in range(5):
        lib.sim_slow(1.0)
    lib.sim_set_faults(1, 0, 0)
    for _ in range(4):
        lib.sim_slow(1.0)
    o = get(lib)
    assert o[8] != 0 and o[9] < 0.5 + 1
    for _ in range(6):
        lib.sim_slow(1.0)
    assert get(lib)[9] == 0 and get(lib)[10] == 0


def test_hardware_trip_input_opens_everything(lib):
    lib.sim_init()
    lib.sim_set_sensors(22.0, 22.0, 6.5415)
    for _ in range(8):
        lib.sim_slow(1.0)
    assert get(lib)[9] == 1
    lib.sim_set_faults(0, 0, 1)
    for _ in range(4):
        lib.sim_slow(1.0)
    o = get(lib)
    assert o[9] == 0 and o[10] == 0
