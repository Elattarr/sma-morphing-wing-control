"""C port vs Python: log closed-loop Python runs, replay the same inputs through the C code, compare sample by sample."""
import os, re, shutil, subprocess, sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
os.environ.update(SMA_D_MM="0.5", SMA_IMAX="3", SMA_SPRING_K="1.25")
import export_replay as ex

pytestmark = pytest.mark.skipif(shutil.which("gcc") is None, reason="needs gcc")
SCENARIOS = ["nominal", "wrong_stress_300", "noisy_sensors", "fault_tc_offset", "fault_stuck_on", "fault_position_frozen", "fault_controller_hang"]


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("cport")
    flags = ["-std=c99", "-Wall", "-Wextra", "-pedantic", "-Werror", "-O2"]
    for name, extra in (("replay", []), ("replay_f", ["-DSMA_FLOAT"])):
        subprocess.run(["gcc", *flags, *extra, str(ROOT / "c/sma.c"), str(ROOT / "c/replay.c"), "-o", str(d / name), "-lm"], check=True)
    logs = {}
    for n, kw in ex.scenarios().items():
        if n in SCENARIOS:
            logs[n] = ex.export(n, d, **kw)
    return d, logs


def run(d, exe, name):
    out = subprocess.run([str(d / exe), str(d / f"{name}.txt")], check=True, capture_output=True, text=True).stdout.split()
    return {out[i]: float(out[i + 1]) for i in range(0, len(out) - 1, 2)}


@pytest.mark.parametrize("name", SCENARIOS)
def test_double_precision_matches_python(built, name):
    d, logs = built
    r = run(d, "replay", name)
    path, py = logs[name]
    assert r["max_dI_ctl"] < 1e-9 and r["max_dI_perm"] < 1e-9 and r["relay_mismatch"] == 0 and r["first_bad"] == -1
    assert int(r["trip"]) == py["trip"]
    cal = py["cal"]
    for k, key in (("as", "As"), ("af", "Af"), ("tfs", "Tfs"), ("tfe", "Tfe")):
        assert r[k] == pytest.approx(getattr(cal, key), abs=1e-6)
    assert r["eps_ext"] == pytest.approx(cal.eps_ext, abs=1e-9)


@pytest.mark.parametrize("name", SCENARIOS)
def test_single_precision_same_decisions(built, name):
    """float: one-sample shifts of threshold decisions are expected; calibration and trip reason must agree"""
    d, logs = built
    r = run(d, "replay_f", name)
    path, py = logs[name]
    assert int(r["trip"]) == py["trip"] and int(r["cal"]) == 1
    cal = py["cal"]
    for k, key in (("as", "As"), ("af", "Af"), ("tfs", "Tfs"), ("tfe", "Tfe")):
        tol = 0.05
        if k == "tfe":                                  # the 99.5 % crossing sits on the slow extension tail: float rounding picks a neighbouring sample
            tol = 1.5 if name == "noisy_sensors" else 0.15
        assert r[k] == pytest.approx(getattr(cal, key), abs=tol)
