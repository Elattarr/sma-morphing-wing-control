"""Identification on simulated bench data (verification only: same model structure as the plant)."""
import os, sys
os.environ.update(SMA_D_MM="0.5", SMA_IMAX="3", SMA_SPRING_K="1.25")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from dataclasses import replace
import pytest
import identify_bench as ib


@pytest.mark.parametrize("name,p,lag", [("nominal", ib.NOM, 0.5), ("slow cooling", replace(ib.NOM, h=5.0), 0.5),
                                        ("high latent", replace(ib.NOM, latent_J_per_mm3=2 * ib.NOM.latent_J_per_mm3), 1.0)])
def test_identification_recovers_parameters(name, p, lag):
    n, r = ib.run_case(name, p, lag)
    assert r is not None
    for k, tol in (("mc", 3.0), ("hS", 3.0), ("RM", 1.0), ("RA", 1.0), ("mL", 10.0)):
        f, t = r[k]
        assert abs(100 * (f - t) / t) < tol, (k, f, t)
