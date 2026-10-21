"""
Standalone SMA wire: thermo-mechanical plant + position PID.

Plant (1-D, lumped, quasi-static mechanics):
  Electrical/thermal : m*c*dT/dt = I^2 R(xi) - h*As*(T - Tinf) + m*dH*dxi/dt
  Phase (hysteresis) : stress-shifted cosine major loops (Liang-Rogers shape) with
                       min/max "play" memory -> nested minor loops, rate independent
                       heating : xi = min(xi, xi_A(T))      (As^s, Af^s)
                       cooling : xi = max(xi, xi_M(T))      (Ms^s, Mf^s)
                       X^s = X + sigma/beta                 (Clausius-Clapeyron)
  Mechanics          : sigma = E(xi) (eps - epsL*xi),  wire in series with bias load
                       sigma*A = F0 + kb*L*eps   ->  eps solved in closed form
  Sensor             : contraction stroke y = L*(eps_cold - eps) [mm] (0 = cold/martensite,
                       positive = wire shortened), noise + quantisation
  Actuator           : current source 0..Imax (PWM-averaged)

Parameters marked [THESIS] come from the master-thesis ANSYS/MEFF work,
[ASSUMED] are placeholders that MUST be replaced by DSC / bench data.
"""
from dataclasses import dataclass, field
import numpy as np


@dataclass
class SMAParams:
    # --- geometry [THESIS: 2 mm wire, 200 mm] -------------------------------
    d: float = 2.0e-3            # wire diameter [m]
    L: float = 0.20              # active length [m]
    # --- material ------------------------------------------------------------
    EA: float = 70e9             # [ASSUMED] austenite modulus [Pa]
    EM: float = 46e9             # [THESIS] martensite modulus [Pa]
    eps_L: float = 0.04          # [THESIS] 3.27 % = sqrt(2/3)*eps_L  ->  eps_L = 4 %
    beta: float = 5.6e6          # [THESIS] Clausius-Clapeyron [Pa/K]
    As: float = 50.8             # [THESIS/MEFF] T0+(R-h*eps_L)/beta
    Af: float = 51.8             # [THESIS/MEFF] T0+R/beta
    Ms: float = 24.0             # [ASSUMED]
    Mf: float = 6.0              # [ASSUMED]
    rho: float = 6450.0          # [ASSUMED] density [kg/m3]
    c: float = 450.0             # [ASSUMED] specific heat [J/kg/K]
    dH: float = 20e3             # [ASSUMED] latent heat [J/kg]
    rho_eA: float = 0.82e-6      # [ASSUMED] resistivity austenite [ohm m]
    rho_eM: float = 0.76e-6      # [ASSUMED] resistivity martensite [ohm m]
    # --- environment / load --------------------------------------------------
    h: float = 30.0              # [ASSUMED] convection coeff [W/m2/K] (light forced air)
    Tinf: float = 22.0           # [THESIS T_load] ambient [degC]
    sigma0: float = 111e6        # bias stress from dead weight / pre-load [Pa]
    kb_L_over_A: float = 190e6   # bias spring stiffness expressed as kb*L/A [Pa]
    # --- electronics ---------------------------------------------------------
    Imax: float = 12.0           # current cap [A]: steady-state T ~100 degC (no T sensor -> cap is the protection)
    noise_mm: float = 0.02       # position sensor noise (1 sigma) [mm]
    quant_mm: float = 0.0        # sensor resolution [mm], 0 = off
    # --- derived -------------------------------------------------------------
    A: float = field(init=False)
    mass: float = field(init=False)
    As_surf: float = field(init=False)

    def __post_init__(self):
        self.A = np.pi * self.d ** 2 / 4
        self.mass = self.rho * self.A * self.L
        self.As_surf = np.pi * self.d * self.L

    @property
    def F0(self):
        return self.sigma0 * self.A

    @property
    def kb(self):
        return self.kb_L_over_A * self.A / self.L


def _cos_curve(T, T0, T1):
    """1 at T<=T0, 0 at T>=T1, half-cosine in between (T0<T1)."""
    if T <= T0:
        return 1.0
    if T >= T1:
        return 0.0
    return 0.5 * (1 + np.cos(np.pi * (T - T0) / (T1 - T0)))


class SMAWire:
    """Plant. step(I, dt) advances one sample and returns the true state."""

    def __init__(self, p: SMAParams, T0=None, xi0=1.0, seed=0):
        self.p = p
        self.T = p.Tinf if T0 is None else T0
        self.xi = xi0
        self.sigma = p.sigma0
        self.rng = np.random.default_rng(seed)
        self.eps = self._strain(self.xi)
        self.eps_ref = self._strain(1.0)       # fully martensitic = stroke zero

    # xi_A / xi_M major-loop curves with stress-shifted temperatures
    def _xiA(self, T):
        s = self.sigma / self.p.beta
        return _cos_curve(T, self.p.As + s, self.p.Af + s)

    def _xiM(self, T):
        s = self.sigma / self.p.beta
        return self._xiM_direct(T, self.p.Mf + s, self.p.Ms + s)   # 0 at Ms^s, 1 at Mf^s

    @staticmethod
    def _xiM_direct(T, Mf_s, Ms_s):
        if T >= Ms_s:
            return 0.0
        if T <= Mf_s:
            return 1.0
        return 0.5 * (1 + np.cos(np.pi * (T - Mf_s) / (Ms_s - Mf_s)))

    def _E(self, xi):
        return self.p.EA + xi * (self.p.EM - self.p.EA)

    def _strain(self, xi):
        p = self.p
        E = self._E(xi)
        k = p.kb_L_over_A
        return (p.sigma0 + E * p.eps_L * xi) / (E - k)

    def _resistance(self, xi):
        p = self.p
        rho_e = p.rho_eA + xi * (p.rho_eM - p.rho_eA)
        return rho_e * p.L / p.A

    def step(self, I, dt, T_force=None):
        """T_force: impose the wire temperature (isothermal sweeps, no latent heat)."""
        p = self.p
        I = float(np.clip(I, 0.0, p.Imax))
        R = self._resistance(self.xi)
        P = I * I * R
        Q = p.h * p.As_surf * (self.T - p.Tinf)
        mc = p.mass * p.c
        T_free = self.T + dt * (P - Q) / mc
        kL = p.dH / p.c                       # adiabatic temperature swing per unit xi
        if T_force is not None:
            T_free, kL = T_force, 0.0
        xi_old = self.xi
        if T_free > self.T:                   # heating branch (xi can only fall)
            if self._xiA(T_free) >= xi_old:
                xi_new = xi_old
            else:
                lo, hi = 0.0, xi_old
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    g = mid - self._xiA(T_free + kL * (mid - xi_old))
                    lo, hi = (mid, hi) if g < 0 else (lo, mid)
                xi_new = 0.5 * (lo + hi)
        elif T_free < self.T:                 # cooling branch (xi can only rise)
            if self._xiM(T_free) <= xi_old:
                xi_new = xi_old
            else:
                lo, hi = xi_old, 1.0
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    g = mid - self._xiM(T_free + kL * (mid - xi_old))
                    lo, hi = (mid, hi) if g < 0 else (lo, mid)
                xi_new = 0.5 * (lo + hi)
        else:
            xi_new = xi_old
        self.T = T_free + kL * (xi_new - xi_old)
        self.xi = xi_new
        self.eps = self._strain(self.xi)
        self.sigma = (p.F0 + p.kb * p.L * self.eps) / p.A      # for next step's shifts
        return self.T, self.xi, self.stroke_mm()

    def stroke_mm(self):
        return 1e3 * self.p.L * (self.eps_ref - self.eps)

    def measure(self):
        p = self.p
        y = self.stroke_mm() + p.noise_mm * self.rng.standard_normal()
        if p.quant_mm > 0:
            y = np.round(y / p.quant_mm) * p.quant_mm
        return y


@dataclass
class PIDGains:
    Kp: float = 3.0      # [A/mm]
    Ki: float = 0.8      # [A/mm/s]
    Kd: float = 0.0      # [A s/mm]
    N: float = 10.0      # derivative filter [1/s]
    Kt: float = 4.0      # anti-windup back-calculation gain [1/s]


class PID:
    """Unipolar (0..Imax) PID, derivative on measurement, back-calculation anti-windup."""

    def __init__(self, g: PIDGains, Imax):
        self.g, self.Imax = g, Imax
        self.i = 0.0
        self.yf = None
        self.d = 0.0

    def reset(self, y):
        self.i, self.yf, self.d = 0.0, y, 0.0

    def update(self, r, y, dt, ff=0.0):
        g = self.g
        if self.yf is None:
            self.yf = y
        # filtered derivative of measurement
        a = g.N * dt / (1 + g.N * dt)
        yf_new = self.yf + a * (y - self.yf)
        self.d = -g.Kd * (yf_new - self.yf) / dt
        self.yf = yf_new
        e = r - y
        u = g.Kp * e + self.i + self.d + ff
        u_sat = float(np.clip(u, 0.0, self.Imax))
        self.i += dt * (g.Ki * e + g.Kt * (u_sat - u))
        return u_sat


def simulate(p: SMAParams, g: PIDGains, ref_fn, t_end, dt=0.01, seed=0, ff_fn=None):
    plant = SMAWire(p, seed=seed)
    pid = PID(g, p.Imax)
    pid.reset(plant.measure())
    n = int(t_end / dt)
    out = {k: np.zeros(n) for k in ("t", "r", "y", "ytrue", "I", "T", "xi")}
    I = 0.0
    for k in range(n):
        t = k * dt
        r = ref_fn(t)
        y_meas = plant.measure()
        I = pid.update(r, y_meas, dt, ff=0.0 if ff_fn is None else ff_fn(r))
        T, xi, ytrue = plant.step(I, dt)
        out["t"][k], out["r"][k], out["y"][k] = t, r, y_meas
        out["ytrue"][k], out["I"][k], out["T"][k], out["xi"][k] = ytrue, I, T, xi
    return out


def step_metrics(t, y, r_final, t_step, y0, band=0.05):
    """rise (10-90 %), overshoot [%], settling (+-band of step size) after t_step."""
    m = t >= t_step
    tt, yy = t[m] - t_step, y[m]
    span = r_final - y0
    if abs(span) < 1e-9:
        return {}
    yn = (yy - y0) / span
    def first(cond):
        idx = np.argmax(cond)
        return tt[idx] if cond.any() else np.nan
    t10, t90 = first(yn >= 0.1), first(yn >= 0.9)
    os_ = max(0.0, (yn.max() - 1) * 100)
    outside = np.where(np.abs(yn - 1) > band)[0]
    ts = tt[outside[-1]] if len(outside) else 0.0
    return {"rise_s": t90 - t10, "overshoot_%": os_, "settle_s": ts}


def static_map(p: SMAParams, n=40, t_settle=400.0, dt=0.05):
    """Steady-state current -> stroke on the HEATING branch (used as feedforward).
    Returns arrays (stroke_mm, I_amps), monotonic in stroke."""
    Is = np.linspace(4.0, p.Imax, n)
    out = []
    for I in Is:
        w = SMAWire(p)
        for _ in range(int(t_settle / dt)):
            w.step(I, dt)
        out.append(w.stroke_mm())
    s, I_ = np.array(out), Is
    # keep the rising part only, make strictly increasing for interpolation
    keep = np.r_[True, np.diff(s) > 1e-6]
    return s[keep], I_[keep]


def make_ff(p: SMAParams, scale=1.0):
    s, I_ = static_map(p)
    return lambda r: 0.0 if r <= 0.05 else scale * float(np.interp(r, s, I_))
