"""heat-only (one-way) controller on the thesis wire: runs the evaluation script in a subprocess (it sets the plant by environment variables)"""
import os, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _run(section):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    out = subprocess.run([sys.executable, str(ROOT / "scripts" / "evaluate_heatonly.py"), section], capture_output=True, text=True, env=env, timeout=600)
    assert out.returncode == 0, out.stderr[-500:]
    return out.stdout


def test_steps_reach_target_without_large_overshoot_and_never_exceed_window():
    txt = _run("steps")
    rows = [l.split() for l in txt.splitlines() if l.strip().startswith("0.00")]
    assert len(rows) == 4
    for r in rows:
        os_pct, sse, tmax = float(r[5]), float(r[6]), float(r[7])
        assert os_pct < 3.0 and abs(sse) < 2.0 and tmax < 60.0
    assert "no trips" in txt


def test_nominal_staircase_reaches_all_targets_with_small_error():
    txt = _run("robust")
    line = next(l for l in txt.splitlines() if l.startswith("nominal"))
    assert "4 of 4" in line
    assert float(re.findall(r"4 of 4\s+([0-9.]+)", line)[0]) < 1.5
