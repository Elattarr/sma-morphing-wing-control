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
from collections import deque
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
    T_abs_max: float = 125.0     # C  absolute temperature limit of the wire/insulation: [ASSUMED] set from the wire datasheet; never exceeded
    T_limit_expand_rate: float = 1.0   # K/s the model-based limit may grow while the wire is still short of its target at the limit
    learn_margin: float = 3.0    # K  after a successful contraction the limit is set to (peak temperature + this), never below the model value
    hold_step: float = 2.0       # K  hold temperature raised (or lowered) by this much when the held position drifts
    hold_drift: float = 0.01     # drift in r (filtered) that triggers a hold-temperature step
    hold_below_As: float = 10.0  # K  the hold temperature never goes above  A_s - this
    overshoot_exit: bool = True  # leave up/down when the position overshoots the target by more than db_start (hold then corrects it)
    sigma_unc: float = 30.0      # MPa, assumed uncertainty of the bias stress: widens the over-temperature limit by 0.146 K/MPa


@dataclass
class Calibration:
    """result of the homing cycle (measured, not modelled)"""
    eps_ext: float        # strain of the cold, fully extended wire
    eps_con: float        # strain of the fully contracted wire
    As: float             # C  temperature at which 0.5 % of the contraction is done (heating, lag corrected)
    Af: float             # C  temperature at which 99.5 % is done
    Tfs: float            # C  temperature at which re-extension starts (cooling): the highest forward start
    Tfe: float            # C  temperature at which re-extension is 99.5 % done
    T_peak: float         # C  highest temperature of the cycle
    t_heat: float         # s
    t_cool: float         # s


class HomingRoutine:
    """Start-up cycle on the cold, extended wire: measure the displacement end points and the four window temperatures.

    wait : thermocouple steady (within 0.3 K over t_wait seconds) and below T_cold_max, displacement not moving: the wire is at
           ambient -> ambient temperature, displacement noise (sets the motion threshold: 5 sigma of a 3 s difference), and the
           displacement of the extended wire
    heat : power ramp  P = Q + m c ramp + I_rate  with a slow integral trim on the heating RATE (not on the temperature): robust to
           thermocouple lag and to a wrong h (time constant tau_rate = 10 s, longer than any sensor lag; outside the window the wire
           heats at `ramp` K/s, inside it the latent heat is absorbed and the window is passed slowly) until the displacement
           has stopped changing for confirm_s seconds after a contraction of at least min_stroke_frac of the expected stroke, or the
           absolute temperature limit is reached -> contracted displacement, As and Af from the 0.5 % and 99.5 % points
    cool : P = 0 until the displacement stops changing (and the wire is near ambient, or has been still for cool_still_max_s)
           -> extended displacement after the cycle, forward start Tfs and end Tfe from the 0.5 % and 99.5 % points
    Motion = displacement change over motion_s seconds above thr_mm (filtered displacement, tau_d). Faults: 'not_cold', 'no_motion_on_heating',
    'no_extension' (the spring cannot return the wire at this ambient temperature).
    """

    def __init__(self, ctl: "CycleController", ramp: float = 1.5, t_wait: float = 5.0, motion_s: float = 2.0, thr_mm: float = 0.05,
                 confirm_s: float = 3.0, min_stroke_frac: float = 0.3, cool_still_s: float = 20.0, cool_still_max_s: float = 90.0,
                 cool_T_margin: float = 4.0, preload_margin_K: float = 2.0, cold_no_ext_s: float = 60.0, t_cool_max: float = 8000.0, t_wait_max: float = 600.0, tau_d: float = 0.5, tau_fast: float = 0.15,
                 T_cold_max: float = 45.0, tau_rate: float = 10.0):
        self.ctl, self.c, self.p = ctl, ctl.c, ctl.p
        self.ramp, self.t_wait, self.thr, self.confirm_s = ramp, t_wait, thr_mm, confirm_s
        self.cool_still_s, self.cool_still_max_s, self.cool_T_margin = cool_still_s, cool_still_max_s, cool_T_margin
        self.preload_margin_K, self.cold_no_ext_s, self._cold = preload_margin_K, cold_no_ext_s, 0.0
        self.t_cool_max, self.t_wait_max, self.tau_d, self.tau_fast = t_cool_max, t_wait_max, tau_d, tau_fast
        self.T_cold_max, self.tau_rate = T_cold_max, tau_rate
        self._Tf, self._Tf_hist, self._I_rate = None, deque(maxlen=int(2.0 / ctl.c.Ts) + 1), 0.0
        self._wait_T, self._wait_d = deque(maxlen=int(t_wait / ctl.c.Ts)), deque(maxlen=int(t_wait / ctl.c.Ts))
        self.min_stroke_mm = min_stroke_frac * ctl.st.u.eps_tr * ctl.p.L * 1e3
        self.phase, self.t, self.done, self.fault, self.cal = "wait", 0.0, False, None, None
        self._df, self._dff, self._hist = None, None, deque(maxlen=int(motion_s / ctl.c.Ts) + 1)
        self._calm = self._still = 0.0
        self._heat, self._cool = [], []
        self._T_peak = -1e9

    def _moving(self, d: float) -> bool:
        a = self.c.Ts / self.tau_d
        self._df = d if self._df is None else self._df + a * (d - self._df)
        af = self.c.Ts / self.tau_fast                      # fast filter: only for the window temperatures of the heating ramp
        self._dff = d if self._dff is None else self._dff + af * (d - self._dff)
        self._hist.append(self._df)
        return len(self._hist) == self._hist.maxlen and abs(self._df - self._hist[0]) > self.thr

    def step(self, T_meas: float, d_meas: float) -> float:
        c, p, ctl = self.c, self.p, self.ctl
        self.t += c.Ts
        moving = self._moving(d_meas)
        self._T_peak = max(self._T_peak, T_meas)
        R = p.resistance(ctl.r_from_disp(d_meas))
        P_max = p.Imax ** 2 * R
        Q = p.h * p.surface * max(T_meas - c.T_amb_est, 0.0)
        if self.phase == "wait":
            self._wait_T.append(T_meas)
            self._wait_d.append(d_meas)
            full = len(self._wait_T) == self._wait_T.maxlen
            steady = False
            if full:
                Tw, dw = np.array(self._wait_T), np.array(self._wait_d)
                h = len(Tw) // 2
                se_d = max(float(np.std(dw)), 1e-6) * np.sqrt(2.0 / h)         # standard error of a difference of half means
                steady = (abs(Tw[h:].mean() - Tw[:h].mean()) < 0.3 and T_meas < self.T_cold_max
                          and abs(dw[h:].mean() - dw[:h].mean()) < max(0.02, 4.0 * se_d))
            if steady:
                c.T_amb_est = float(np.mean(Tw))                                  # the wire is at ambient
                s_f = float(np.std(dw)) * np.sqrt(c.Ts / (2.0 * self.tau_d))      # std of the filtered displacement
                self.thr = max(self.thr, 5.0 * np.sqrt(2.0) * s_f)                 # motion = change over motion_s above 5 sigma
                self._d_e0, self._goal, ctl.Tcmd, self.phase = self._df, T_meas, None, "heat"
            elif self.t > self.t_wait_max:
                self.fault, self.done = "not_cold", True
            return 0.0
        if self.phase == "heat":
            self._heat.append((self.t, T_meas, self._dff))
            stroke = abs(self._df - self._d_e0)
            self._still = 0.0 if moving else self._still + c.Ts
            if T_meas >= c.T_abs_max - 3.0 or (stroke > self.min_stroke_mm and self._still >= self.confirm_s):
                if stroke < self.min_stroke_mm:
                    self.fault, self.done = "no_motion_on_heating", True
                    return 0.0
                self._d_c, self.phase, self._still, self._t_cool0 = self._df, "cool", 0.0, self.t
                return 0.0
            self._Tf = T_meas if self._Tf is None else self._Tf + c.Ts / 1.0 * (T_meas - self._Tf)      # smoothed temperature (tau 1 s)
            self._Tf_hist.append(self._Tf)
            if len(self._Tf_hist) == self._Tf_hist.maxlen:
                rate = (self._Tf_hist[-1] - self._Tf_hist[0]) / (2.0)                                # K/s over 2 s
                self._I_rate += c.Ts * p.heat_capacity / self.tau_rate * (self.ramp - rate)
            P_un = Q + p.heat_capacity * self.ramp + self._I_rate
            P = float(np.clip(P_un, 0.0, P_max))
            if P != P_un:                                                                               # anti-windup on the trim
                self._I_rate += P - P_un
            return float(np.sqrt(P / R)) if T_meas <= c.T_abs_max else 0.0
        # cool
        self._cool.append((self.t, T_meas, self._df))
        started = abs(self._df - self._d_c) > self.thr
        n20 = int(self.cool_still_s / c.Ts)
        drifting = len(self._cool) > n20 and abs(self._df - self._cool[-1 - n20][2]) > self.thr    # slow tail: below thr per 2 s, not per 20 s
        self._still = 0.0 if (moving or drifting or not started) else self._still + c.Ts     # stillness counts only after the extension began
        self._cold = self._cold + c.Ts if (not started and T_meas - c.T_amb_est < self.cool_T_margin) else 0.0
        settled = len(self._cool) > n20 and abs(self._cool[-1][1] - self._cool[-1 - n20][1]) < 0.2     # < 0.01 K/s over 20 s
        if started and self._still >= self.cool_still_s and (T_meas - c.T_amb_est < self.cool_T_margin or
                                                              (self._still >= self.cool_still_max_s and settled)):
            self._finish()
        elif self._cold >= self.cold_no_ext_s or self.t - self._t_cool0 > self.t_cool_max:        # at ambient for a minute and the wire has not moved
            self.fault, self.done = "no_extension", True
        return 0.0

    def _finish(self):
        L, n_persist = self.p.L * 1e3, max(int(0.5 / self.c.Ts), 1)
        H, C = np.array(self._heat), np.array(self._cool)
        fh = (H[:, 2] - self._d_e0) / (self._d_c - self._d_e0)
        d_e2 = C[-1, 2]
        fc = (C[:, 2] - self._d_c) / (d_e2 - self._d_c)

        def at(A, f, x, delay=0.0):
            """temperature at the first time the fraction f is >= x for 0.5 s in a row, minus the filter delay (looked up in time)"""
            run = np.convolve((f >= x).astype(int), np.ones(n_persist, int), "valid") == n_persist
            i = int(np.argmax(run)) if run.any() else len(f) - 1
            return float(np.interp(A[i, 0] - delay, A[:, 0], A[:, 1]))

        # crossing fractions: 0.5 % / 99.5 % of the stroke, but never below the measured noise floor (the motion threshold, 5 sqrt(2) sigma of
        # the filtered displacement): with a noisy position sensor a fixed 0.5 % is crossed by noise long before the real motion starts
        xh = min(max(0.005, self.thr / max(abs(self._d_c - self._d_e0), 1e-9)), 0.25)
        xc = min(max(0.005, self.thr / max(abs(d_e2 - self._d_c), 1e-9)), 0.25)
        # first-order filter delay: tau_fast at the start of a ramp-like motion, tau_fast ln(1/(1-x)) at the end of a step-like one
        self.cal = Calibration(eps_ext=d_e2 / L, eps_con=self._d_c / L,
                               As=at(H, fh, xh, self.tau_fast), Af=at(H, fh, 1.0 - xh, self.tau_fast * np.log(1.0 / xh)),
                               Tfs=at(C, fc, xc), Tfe=at(C, fc, 1.0 - xc),
                               T_peak=self._T_peak, t_heat=len(H) * self.c.Ts, t_cool=len(C) * self.c.Ts)
        self.done = True
        if self.cal.Tfe < self.c.T_amb_est + self.preload_margin_K:          # the extension ended at ambient: the wire did not return fully
            self.fault = "preload_too_low"
            self.cal = None


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
        self.T_limit_model = self.T_limit                                   # the limit can grow in a move and is learned back after it
        self.T_lim = self.T_limit
        self.T_hold = self.T_safe
        self.T_hold_max = max(self.T_safe, self.As - self.c.hold_below_As)
        self.T_peak = -1e9
        self.r_hold_ref = 0.0
        self.cal = None
        self.T_motion = self.Tfs + np.sqrt(2.0 / 3.0) / M.beta * self.c.sigma_unc + 2.0   # highest T where motion may start the down loop

    def apply_calibration(self, cal: "Calibration"):
        """replace the windows and the displacement end points by the values measured in the homing cycle"""
        c, p = self.c, self.p
        self.cal = cal
        self.As, self.Af, self.Tfs, self.Tfe = cal.As, cal.Af, cal.Tfs, cal.Tfe
        self.W = max(self.Af - self.As, 0.5)
        self.m_eff = p.latent_total + p.heat_capacity * self.W
        mid = 0.5 * (self.Tfs + self.Tfe)
        self.rate_dn = c.dn_frac * p.h * p.surface * max(mid - c.T_amb_est, 1.0) / self.m_eff
        self.T_pre = self.As - c.margin_up
        self.T_safe = self.Tfs + c.T_safe_above_fs
        self.T_limit_model = self.T_lim = self.T_limit = min(self.Af + c.T_limit_above_Af, c.T_abs_max)
        self.T_hold = self.T_safe
        self.T_hold_max = max(self.T_safe, self.As - c.hold_below_As)
        self.T_motion = self.Tfs + 2.0
        self._enter("hold")

    # ---------------------------------------------------------------- position from displacement
    def r_from_disp(self, d_mm: float) -> float:
        eps = d_mm / (self.p.L * 1e3)
        cal = self.cal
        if cal is not None:                      # two-point correction: measured end points -> the model's end points, model shape in between
            m1, m0 = self.st.strain(1.0), self.st.strain(0.0)
            eps = m1 + (eps - cal.eps_ext) * (m0 - m1) / (cal.eps_con - cal.eps_ext)
        return 1.0 - self.st.xi_from_strain(eps)

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
        self.r_hold_ref, self.T_peak = self.r_f, -1e9
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
        self.T_peak = max(self.T_peak, T_meas)
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
                if self.in_loop and T_meas > self.T_lim - 1.0 and r_target - self.r_f > c.db_start:
                    # still short of the target at the limit: the real window is wider or higher than the model says
                    self.T_lim = min(self.T_lim + c.T_limit_expand_rate * c.Ts, c.T_abs_max)
                if T_meas > self.T_lim:
                    P = 0.0
                if abs(self.rref - r_target) < 1e-9 and abs(r_target - r) < tol:
                    self.T_lim = max(self.T_limit_model, min(self.T_peak + c.learn_margin, c.T_abs_max))
                    # a part of the wire that needed d K more than A_f to recover also starts its forward transformation d K higher:
                    # hold at least d above the model's safe temperature (never above T_hold_max)
                    self.T_hold = min(max(self.T_hold, self.T_safe + max(self.T_peak - self.Af, 0.0)), self.T_hold_max)
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
            d = self.r_f - self.r_hold_ref
            if d < -c.hold_drift:                                     # re-extending: some of the wire still transforms at this temperature
                self.T_hold, self.r_hold_ref = min(self.T_hold + c.hold_step, self.T_hold_max), self.r_f
            elif d > c.hold_drift:                                    # creeping back: lower again, never below the model value
                self.T_hold, self.r_hold_ref = max(self.T_hold - c.hold_step, self.T_safe), self.r_f
            P = self._temp_loop(self.T_hold, T_meas, Q, P_max, c.Kp_hold, c.Ki_hold)
        if T_meas > self.T_lim or T_meas > c.T_abs_max:                    # over-temperature cut-off in every mode
            P = 0.0
        return float(np.sqrt(P / R))
