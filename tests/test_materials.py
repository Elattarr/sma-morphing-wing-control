"""Other materials, unchanged gains: the controller (with homing) must reach all four targets for different MEFF constants."""
import os, sys
os.environ.update(SMA_D_MM="0.5", SMA_IMAX="3", SMA_SPRING_K="1.25")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from dataclasses import replace
import pytest
import evaluate_materials as em
ev = em.ev


@pytest.mark.parametrize("label,kw", [("wide window h=1400", dict(h=1400.0)), ("beta 4.5", dict(beta=4.5)), ("R 220", dict(R=220.0))])
@pytest.mark.parametrize("own", [False, True])
def test_other_alloy_reaches_all_targets(label, kw, own):
    M = replace(em.M0, **kw)
    p_true = replace(ev.NOM, meff=M)
    sig = em.preload(M)
    L, seg, home = ev.sequence(p_true=p_true, sigma_true=sig, sigma_ctl=sig, p_ctl=p_true if own else None, label=label)
    assert home["fault"] is None and len(seg) == 4
    assert not home["guard"].mon.tripped
