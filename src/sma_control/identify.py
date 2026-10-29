"""Parameter identification from bench data (current I, voltage V, thermocouple T_meas, displacement d) of ONE wire.

The controller needs a handful of numbers that datasheets give only loosely; three short tests measure them:

  Test A  THERMAL STEP (wire extended, cold, stays below A_s): constant current, cut when T_meas - T_amb reaches a limit, then free cooling.
          Energy balance with a first-order thermocouple lag:
              m c dT/dt = I^2 R - h S (T - T_amb)          T_meas = lag(T, tau_lag)
          fitted by nonlinear least squares on the whole record  ->  m c, h S, tau_lag.
          R = mean(V / I) over the heating part  ->  R_M (martensite resistance).
  Test B  LATENT HEAT (the homing heating record, wire passes through the whole transformation): with m c, h S, tau_lag known,
          the temperature response is LINEAR in the latent heat m Lambda, so it is a one-parameter linear least squares:
              m c dT/dt = I^2 R(r) - h S (T - T_amb) - m Lambda dr/dt          r(t) from the displacement
  Test C  RESISTANCE of the contracted wire: mean(V/I) where the displacement says the wire is fully recovered  ->  R_A.

What is NOT identified: the transformation windows and the strains (homing measures them, see HomingRoutine), and the stress (preload;
measured with a load cell or from the spring). Everything here is verified on simulated data only.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares


@dataclass
class ThermalFit:
    mc: float            # J/K
    hS: float            # W/K
    tau_lag: float       # s
    R: float             # ohm  (resistance during the test)
    rms: float           # K    residual rms
    T_amb: float


def _simulate(theta, t, P, T_amb, T0):
    """first-order wire + first-order thermocouple, exact discretisation of the lag, Euler for the wire (dt << tau)"""
    mc, hS, tau = theta
    n = len(t)
    T, Tm = np.empty(n), np.empty(n)
    T[0] = Tm[0] = T0
    for k in range(1, n):
        dt = t[k] - t[k - 1]
        T[k] = T[k - 1] + dt * (P[k - 1] - hS * (T[k - 1] - T_amb)) / mc
        a = dt / (tau + dt)
        Tm[k] = Tm[k - 1] + a * (T[k] - Tm[k - 1])
    return Tm, T


def fit_thermal(t, I, V, Tm, T_amb=None, guess=None) -> ThermalFit:
    """Test A. t, I, V, Tm: arrays of the whole record (heating + cooling). Returns the fitted parameters."""
    t, I, V, Tm = map(np.asarray, (t, I, V, Tm))
    on = I > 0.5 * I.max()
    R = float(np.mean(V[on] / I[on]))
    T_amb = float(np.mean(Tm[: max(int(8.0 / (t[1] - t[0])), 5)])) if T_amb is None else T_amb          # record starts with >= 8 s at rest
    P = (I * I) * R
    # starting values: initial heating slope -> m c, cooling decay -> h S / m c
    n0 = max(int(1.5 / (t[1] - t[0])), 5)
    j = int(np.argmax(on))
    slope = max((Tm[j + 2 * n0] - Tm[j + n0]) / (t[j + 2 * n0] - t[j + n0]), 1e-3)
    mc0 = float(P[j] / slope) if guess is None else guess[0]
    hS0 = float(P[j] * 0.1 / max(Tm.max() - T_amb, 1.0)) if guess is None else guess[1]

    def res(x):
        th = np.exp(x)
        Tm_pred, _ = _simulate(th, t, P, T_amb, T_amb)
        return Tm_pred - Tm

    x0 = np.log([mc0, hS0, 0.5])
    sol = least_squares(res, x0, bounds=(np.log([1e-3, 1e-5, 0.01]), np.log([1e3, 10.0, 20.0])), x_scale=1.0)
    mc, hS, tau = np.exp(sol.x)
    return ThermalFit(mc=float(mc), hS=float(hS), tau_lag=float(tau), R=R, rms=float(np.sqrt(np.mean(sol.fun ** 2))), T_amb=T_amb)


def fit_latent(t, I, V, Tm, d, th: ThermalFit, d_ext: float, d_con: float, t_end: float | None = None) -> float:
    """Test B. The record is the homing HEATING phase (cold -> wire fully contracted). Returns m*Lambda [J].
    The temperature is linear in m Lambda:  Tm = Tm0 + m Lambda * Tm1, solved in closed form by least squares."""
    t, I, V, Tm, d = map(np.asarray, (t, I, V, Tm, d))
    if t_end is not None:
        k = t <= t_end
        t, I, V, Tm, d = t[k], I[k], V[k], Tm[k], d[k]
    dt = float(t[1] - t[0])
    w = max(int(1.0 / dt), 1)
    r = np.clip((d - d_ext) / (d_con - d_ext), 0.0, 1.0)
    r = np.convolve(np.pad(r, (w // 2, w - w // 2 - 1), mode="edge"), np.ones(w) / w, mode="valid")        # 1 s moving average
    r = np.maximum.accumulate(r)                                                                          # heating: r only grows
    rdot = np.gradient(r, dt)
    Rw = np.where(I > 0.2, V / np.maximum(I, 1e-9), th.R)
    P = I * I * Rw
    # response to the known forcing, and to a unit latent forcing  -rdot  (same linear system)
    Tm0, _ = _simulate((th.mc, th.hS, th.tau_lag), t, P, th.T_amb, th.T_amb)
    Tm1, _ = _simulate((th.mc, th.hS, th.tau_lag), t, -rdot, th.T_amb, th.T_amb)
    Tm1 = Tm1 - th.T_amb                                                                                   # unit forcing starts from zero
    # _simulate returned T_amb + response because T starts at T_amb and the loss term is relative to T_amb: remove the offset
    a = np.sum(Tm1 * (Tm - Tm0)) / np.sum(Tm1 * Tm1)
    return float(a)


def fit_resistance_contracted(V, I, d, d_ext: float, d_con: float) -> float:
    """Test C. R of the fully recovered wire: samples with r > 0.98 and I > 0.3 A"""
    V, I, d = map(np.asarray, (V, I, d))
    r = (d - d_ext) / (d_con - d_ext)
    k = (r > 0.98) & (I > 0.3)
    return float(np.mean(V[k] / I[k])) if k.any() else float("nan")
