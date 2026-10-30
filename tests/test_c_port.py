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


# ---- heat-only controller: exported in its own process (the plant is chosen by environment variables at import) ----
HO = ["ho_nominal", "ho_spread4", "ho_noisy", "ho_short_stroke", "ho_fault_tc_offset", "ho_fault_stuck_on"]


@pytest.fixture(scope="module")
def built_ho(tmp_path_factory):
    d = tmp_path_factory.mktemp("cport_ho")
    flags = ["-std=c99", "-Wall", "-Wextra", "-pedantic", "-Werror", "-O2"]
    for name, extra in (("replay", []), ("replay_f", ["-DSMA_FLOAT"])):
        subprocess.run(["gcc", *flags, *extra, str(ROOT / "c/sma.c"), str(ROOT / "c/replay.c"), "-o", str(d / name), "-lm"], check=True)
    subprocess.run([sys.executable, str(ROOT / "scripts/export_replay_heatonly.py"), str(d)], check=True, capture_output=True,
                   env={k: v for k, v in os.environ.items() if not k.startswith("SMA_")})
    return d


@pytest.mark.parametrize("name", HO)
def test_heat_only_double_precision_identical(built_ho, name):
    r = run(built_ho, "replay", name)
    assert r["max_dI_ctl"] < 1e-9 and r["max_dI_perm"] < 1e-9 and r["relay_mismatch"] == 0 and r["first_bad"] == -1
    assert int(r["ctl_fault"]) == 0


def test_heat_only_learner_path_is_replayed(built_ho):
    r = run(built_ho, "replay", "ho_short_stroke")
    assert int(r["ctl_cal"]) == 1 and int(r["ctl_sat"]) == 1 and 52.0 < r["ctl_af"] < 66.0 and int(r["trip"]) == 0


@pytest.mark.parametrize("name", HO)
def test_heat_only_single_precision_same_decisions(built_ho, name):
    ref = run(built_ho, "replay", name)
    r = run(built_ho, "replay_f", name)
    assert int(r["trip"]) == int(ref["trip"]) and int(r["ctl_cal"]) == int(ref["ctl_cal"]) and int(r["ctl_fault"]) == int(ref["ctl_fault"])
    assert r["max_dI_ctl"] < 0.05
