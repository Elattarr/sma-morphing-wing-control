import numpy as np
from sma_control import SMAParams
from sma_control.spring import SMAWireSpring, BiasSpring
from sma_control.elements import SMAWireElements, element_offsets

P, SP = SMAParams(), BiasSpring(250.0, 20.0)


def run(w, I_of, t_end, dt=0.05):
    log = []
    for k in range(int(t_end / dt)):
        w.step(I_of(k * dt, w), dt, 22.0)
        log.append((w.T, w.xi))
    return np.array(log)


def heat_then_cool(w, T_hot=110.0, t_end=1500.0):
    hot = [True]

    def I_of(t, w):
        if w.T >= T_hot:
            hot[0] = False
        return 15.0 if hot[0] else 0.0
    return run(w, I_of, t_end)


def test_zero_spread_equals_spring_plant():
    a = SMAWireSpring(P, SP, 22.0, 1.0)
    b = SMAWireElements(P, SP, 22.0, 1.0, spread_K=0.0)
    la, lb = heat_then_cool(a), heat_then_cool(b)
    assert np.abs(la - lb).max() < 1e-6
    assert la[-1, 1] > 0.99


def test_offsets_are_symmetric_gaussian_plus_tail():
    o = element_offsets(21, 3.0, 0.0, 0.0)
    assert abs(o.mean()) < 1e-12 and abs(o.std() - 3.0) < 0.35
    o2 = element_offsets(20, 2.0, 0.25, 20.0)
    assert (o2 > 4).sum() == 4 and o2.max() < 20.0 and len(o2) == 20


def test_window_is_broader_by_the_spread():
    w0 = SMAWireSpring(P, SP, 22.0, 1.0)
    w1 = SMAWireElements(P, SP, 22.0, 1.0, spread_K=3.0)
    def width(w):
        L = heat_then_cool(w)
        T, xi = L[:, 0], L[:, 1]
        hot = np.argmax(T)                                             # end of heating
        h = np.where((np.arange(len(T)) <= hot))[0]
        t95 = T[h][np.argmax(xi[h] < 0.95)]
        t05 = T[h][np.argmax(xi[h] < 0.05)]
        return t05 - t95
    assert width(w1) > width(w0) + 6.0                                 # 2 x 1.64 x 3 K = 9.8 K extra, minus coarse quantiles


def test_sub_loop_return_point_memory():
    """heat to the middle of the window, cool, reheat: the wire stays frozen until the previous peak temperature is passed again"""
    w = SMAWireElements(P, SP, 22.0, 1.0, spread_K=3.0)
    T_mid = w.st.T_rec(0.5)
    peak = -1e9
    for _ in range(int(800 / 0.05)):                                   # heat to T_mid, then cool back to 40 C
        w.step(15.0 if (w.T < T_mid and peak < T_mid) else 0.0, 0.05, 22.0)
        peak = max(peak, w.T)
        if peak >= T_mid and w.T < 40.0:
            break
    xi_after = w.xi
    assert 0.2 < xi_after < 0.8
    for _ in range(int(2000 / 0.05)):                                  # reheat slowly (temperature-controlled) to 1 K below the peak
        w.step(min(max(3.0 + 0.8 * (peak - 1.0 - w.T), 0.0), 6.0) if w.T < peak - 1.0 else 0.0, 0.05, 22.0)
        if w.T >= peak - 1.0:
            break
    assert abs(w.xi - xi_after) < 5e-3                                 # frozen below the previous peak
    for _ in range(int(60 / 0.05)):                                    # pass the peak: transformation resumes
        w.step(15.0, 0.05, 22.0)
        if w.T > peak + 3.0:
            break
    assert w.xi < xi_after - 0.03
