"""Bidirectional position controller for the one-way wire under a constant bias stress sigma (cycle_plant.py).

Physics it relies on (MEFF, one-way, constant stress):
  * contraction (r up)    happens only in the recovery window   [A_s(sigma), A_f(sigma)]            (hot, e.g. 87-88 C at 250 MPa)
  * extension  (r down)   happens only in the forward window    [T_fs - W, T_fs]  (T_fs = forward start) (cool, e.g. 27.7-28.7 C)
  * between the two windows xi is frozen (about 60 K dead zone): position is held for free.
  * inside either window the heat balance gives the same integrator in r:
        (m*Lambda + m*c*W) dr/dt = P - h S (T - T_amb)          (P = I^2 R)
    so one position law serves both directions:   P = Q + m*Lambda*(drref/dt + Kr*e + Ki_r*int e),  Q = h S (T - T_amb)
  * extension is limited by heat loss alone:   |dr/dt| <= Q(T_fs) / (m*Lambda)   (about 0.0026 per s here): cannot be commanded faster.

Modes
  hold : temperature loop to T_safe = T_fs + 3 K. Position is frozen there; the wire is kept just above the forward window
         so it cannot re-extend by itself (costs about h S (T_safe - T_amb) = 0.4 W).
  up   : temperature loop to A_s - margin, then the position law with a rate-limited rising reference.
  down : P = 0 (cool) until T < T_fs + margin, then the position law with a falling reference.
  In up and down the position law also starts as soon as the position has moved (the believed window can be wrong
  because sigma is uncertain): preheat and passive cooling are never continued once the wire is already moving.
Position feedback is the displacement; r is recovered from it with the known stress:   d/L = sigma/E(xi) + eps_tr*xi,  r = 1 - xi.
The controller knows sigma and the nominal MEFF numbers: a wrong sigma moves the windows (0.146 K/MPa), see scripts/control_bias_cycle.py.

Bias spring (K_spring > 0, see spring.py): the stress is sigma(xi) = sigma_ext + K (eps_ref - eps) instead of a constant, so
  * r comes from the displacement through the spring law (SpringStress.xi_from_strain),
  * the recovery window is [T_rec(1), T_rec(0)] and the forward window [T_fwd(1), T_fwd(0)], each wider than the free-wire 1 K,
  * inside a window the heat balance has the effective width W_eff = T_rec(0) - T_rec(1):
        (m*Lambda + m*c*W_eff) dr/dt = P - Q
K_spring = 0 gives back the constant-stress controller above (sigma is then sigma_ext).
"""
from dataclasses import dataclass
import numpy as np
from .params import SMAParams
from .loading import UniaxialMEFF
from .spring import SpringStress


@dataclass
class CycleCfg:
    Ts: float = 0.02
    Kp: float = 5.0              # W/K     temperature loop
    Ki: float = 10.0             # W/K/s
    k_soft: float = 3.0          # 1/s
    dTdt_max: float = 20.0       # K/s
    margin_up: float = 1.0       # K below A_s(sigma) where the preheat stops
    margin_dn: float = 0.3       # K above T_fs where the down loop takes over
    T_safe_above_fs: float = 3.0
    Kr: float = 2.0              # 1/s     position loop bandwidth, contraction
    Ki_r: float = 0.3
    Kr_dn: float = 0.05          # 1/s     extension is heat-loss limited (about 100 s time scale): slow loop,
    Ki_r_dn: float = 0.0005      #         and a filtered position, so sensor noise is not rectified by the P >= 0 limit
    tau_filt: float = 3.0        # s       position filter for the extension loop and for the start decisions
    Kp_hold: float = 0.5         # W/K     hold needs no speed: low gains keep sensor noise out of the current
    Ki_hold: float = 0.05        # W/K/s
    drdt_up: float = 0.3         # 1/s
    dn_frac: float = 0.8         # fraction of the physical maximum extension rate used as reference rate
    Kaw: float = 5.0
    db_start: float = 0.02       # position error that starts a move
    db_done: float = 0.01        # position error that ends a move (reference already at the target)
    T_amb_est: float = 22.0
    T_limit_above_Af: float = 3.0
    overshoot_exit: bool = True  # leave up/down when the position overshoots the target by more than db_start (hold then corrects it)
    sigma_unc: float = 30.0      # MPa, assumed uncertainty of the bias stress: widens the over-temperature limit by 0.146 K/MPa


def cfg_for_wire(p: SMAParams, **overrides) -> CycleCfg:
    """CycleCfg with the temperature-loop gains (W/K, W/K/s) scaled by the wire's heat capacity m*c, relative to the 2 mm wire
    they were tuned on: the loop time constant m*c/Kp must stay the same, otherwise a thin wire (0.5 mm: 16 times less mass)
    is driven at a loop bandwidth beyond the sample rate. The contraction position-loop gains are normalised by m*Lambda already."""
    ref = SMAParams()
    k = p.heat_capacity / ref.heat_capacity
    base = CycleCfg()
    scaled = dict(Kp=base.Kp * k, Ki=base.Ki * k, Kp_hold=base.Kp_hold * k, Ki_hold=base.Ki_hold * k)
    # extension is limited by heat loss, so its time scale follows the cooling time constant tau = m c / (h S):
    # the slow extension loop and its filter are scaled with it (tau_ratio = 0.25 for the 0.5 mm wire)
    tr = p.cooling_tau / ref.cooling_tau
    scaled.update(Kr_dn=base.Kr_dn / tr, Ki_r_dn=base.Ki_r_dn / tr ** 2, tau_filt=base.tau_filt * tr)
    scaled.update(overrides)
    return CycleCfg(**scaled)


class CycleController:
    def __init__(self, p: SMAParams, sigma: float, cfg: CycleCfg | None = None, K_spring: float = 0.0):
        """sigma = stress with the wire fully extended (the whole stress if K_spring = 0); K_spring in MPa per unit strain"""
        self.p, self.c, self.sigma = p, cfg or CycleCfg(), sigma
        self.u = UniaxialMEFF(p.meff)
        M, u = p.meff, self.u
        self.st = SpringStress(u, sigma, K_spring)
        self.As, self.Af = self.st.T_rec(1.0), self.st.T_rec(0.0)          # recovery starts at xi = 1, ends at xi = 0
        self.Tfs = self.st.T_fwd(0.0)                                        # highest forward start (wire contracted, highest stress)
        self.Tfe = self.st.T_fwd(1.0)                                        # temperature at which re-extension is complete
        self.W = self.Af - self.As                                           # effective window width (= h eps_L / beta for K = 0)
        self.m_eff = p.latent_total + p.heat_capacity * self.W              # m*Lambda + m*c*W_eff
        self.T_limit = self.Af + self.c.T_limit_above_Af + np.sqrt(2.0 / 3.0) / M.beta * self.c.sigma_unc
        self.T_pre = self.As - self.c.margin_up
        self.T_safe = self.Tfs + self.c.T_safe_above_fs
        mid = 0.5 * (self.Tfs + self.Tfe)
        self.rate_dn = self.c.dn_frac * p.h * p.surface * (mid - self.c.T_amb_est) / self.m_eff
        self.mode = "hold"
        self.in_loop = False
        self.Tcmd = None
        self.I_a = self.I_b = 0.0
        self.rref = 0.0
        self.r_est = 0.0
        self.r_f = None
        self.r_entry = 0.0
        self.rf_entry = 0.0
        self.T_motion = self.Tfs + np.sqrt(2.0 / 3.0) / M.beta * self.c.sigma_unc + 2.0   # highest T where motion may start the down loop

    # ---------------------------------------------------------------- position from displacement
    def r_from_disp(self, d_mm: float) -> float:
        return 1.0 - self.st.xi_from_strain(d_mm / (self.p.L * 1e3))

    # ---------------------------------------------------------------- loops
    def _temp_loop(self, goal, T_meas, Q, P_max, Kp=None, Ki=None):
        c = self.c
        Kp = c.Kp if Kp is None else Kp
        Ki = c.Ki if Ki is None else Ki
        if self.Tcmd is None:
            self.Tcmd = T_meas
        step = np.clip(c.k_soft * (goal - self.Tcmd), -c.dTdt_max, c.dTdt_max) * c.Ts
        if abs(goal - self.Tcmd) > c.dTdt_max / c.k_soft:           # constant-rate ramp far from the goal
            step = np.sign(goal - self.Tcmd) * c.dTdt_max * c.Ts
        self.Tcmd += step
        e = self.Tcmd - T_meas
        P_un = Q + Kp * e + self.I_a
        P = float(np.clip(P_un, 0.0, P_max))
        self.I_a += c.Ts * (Ki * e + c.Kaw * (P - P_un))
        return P

    def _pos_loop(self, r_target, r_meas, Q, P_max, rate_up, Kr, Ki_r):
        c, p = self.c, self.p
        if not self.in_loop:
            self.in_loop, self.rref, self.I_b = True, r_meas, 0.0
        step = float(np.clip(r_target - self.rref, -self.rate_dn * c.Ts, rate_up * c.Ts)) if rate_up else 0.0
        self.rref += step
        e = self.rref - r_meas
        mL = self.m_eff
        P_un = Q + mL * (step / c.Ts + Kr * e + self.I_b)
        P = float(np.clip(P_un, 0.0, P_max))
        self.I_b += c.Ts * (Ki_r * e + c.Kaw * (P - P_un) / mL)
        return P

    def _enter(self, mode):
        self.r_entry = self.r_est
        self.r_f = self.r_est                                   # restart the filter at the current position
        self.rf_entry = self.r_f
        self.mode, self.in_loop, self.Tcmd = mode, False, None
        self.I_a = self.I_b = 0.0

    # ---------------------------------------------------------------- one control step
    def step(self, r_target: float, T_meas: float, d_meas_mm: float) -> float:
        c, p = self.c, self.p
        r = self.r_from_disp(d_meas_mm)
        self.r_est = r
        self.r_f = r if self.r_f is None else self.r_f + c.Ts / c.tau_filt * (r - self.r_f)
        R = p.resistance(r)
        P_max = p.Imax ** 2 * R
        Q = p.h * p.surface * max(T_meas - c.T_amb_est, 0.0)
        if self.mode == "hold":
            if r_target > self.r_f + c.db_start:
                self._enter("up")
            elif r_target < self.r_f - c.db_start:
                self._enter("down")
        P = 0.0
        # at a bound (r = 0 or 1) the estimate cannot cross it: a wrong sigma biases r_est by 0.07 % strain per 25 MPa,
        # so accept 0.03 there; inside the range use the tight tolerance
        tol = c.db_done if 0.03 < r_target < 0.97 else 0.03
        if self.mode == "up":
            if not self.in_loop and T_meas < self.T_pre - 0.1 and r < self.r_entry + c.db_start:
                P = self._temp_loop(self.T_pre, T_meas, Q, P_max)
            else:
                P = self._pos_loop(r_target, r, Q, P_max, c.drdt_up, c.Kr, c.Ki_r)
                if T_meas > self.T_limit:
                    P = 0.0
                if abs(self.rref - r_target) < 1e-9 and abs(r_target - r) < tol:
                    self._enter("hold")
                elif c.overshoot_exit and self.r_f > r_target + c.db_start:          # overshot beyond the tolerance: stop heating, let hold decide (extend)
                    self._enter("hold")
        elif self.mode == "down":
            if not self.in_loop and (T_meas > self.T_motion or (T_meas > self.Tfs + c.margin_dn and self.r_f > self.rf_entry - 0.5 * c.db_start)):
                P = 0.0
            else:
                P = self._pos_loop(r_target, self.r_f, Q, P_max, c.drdt_up, c.Kr_dn, c.Ki_r_dn)
                if abs(self.rref - r_target) < 1e-9 and abs(self.r_f - r_target) < tol:
                    self._enter("hold")
                elif c.overshoot_exit and self.in_loop and self.r_f < r_target - c.db_start:   # extended too far: back to hold (hold decides to contract)
                    self._enter("hold")
        else:                                                         # hold: keep the wire just above the forward window
            P = self._temp_loop(self.T_safe, T_meas, Q, P_max, c.Kp_hold, c.Ki_hold)
        if T_meas > self.T_limit:                                          # over-temperature cut-off in every mode
            P = 0.0
        return float(np.sqrt(P / R))
