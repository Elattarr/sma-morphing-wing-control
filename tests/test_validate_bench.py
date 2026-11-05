"""validate_bench.py must accept a log from the model itself and reject one from a wire with 1.6x the convection (checks the script, not the wire)."""
import subprocess, sys
from pathlib import Path

def test_selftest_accepts_and_rejects():
    r = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / "scripts" / "validate_bench.py"), "--selftest"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    a, b = r.stdout.split("== selftest: true h and c x1.6")
    assert "VALIDATED for this run" in a and "NOT validated" in b
