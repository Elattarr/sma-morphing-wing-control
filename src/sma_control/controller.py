"""Baseline controller for the standalone one-way wire.

Inside the 1 K window the plant is an integrator in recovery r (latent heat dominates):

    m*Lambda * dr/dt = P - h*S*(T - T_amb)          P = I^2 R

Temperature is a poor feedback there (0.1 K noise = 10 % of the stroke), so the controller has two phases:

  A  preheat   : temperature loop (P_ff + Kp*e + Ki*int e) to T_hold = As - margin; rate-limited, soft landing.
                 Position is blind in this phase (no output for the first ~30 K).
  B  position  : feedback-linearised integrator, in the power domain
                     P = P_loss(T_meas) + m*Lambda*( drref/dt + Kr*(rref - r_meas) + Ki_r*int(rref - r_meas) )
                 -> dr/dt = Kr*e, first order, no overshoot by design (recovery cannot be undone).
Current I = sqrt(P / R(r_est)), clipped to [0, Imax]; back-calculation anti-windup on both integrators.
"""
from dataclasses import dataclass
import numpy as np
from .params import SMAParams, RecoveryMap


@dataclass
class CtrlCfg:
    Ts: float = 0.01
    # phase A
    Kp: float = 5.0             # W/K
    Ki: float = 10.0            # W/K/s
    k_soft: float = 3.0         # 1/s  soft landing of the temperature target
    dTdt_max: float = 20.0      # K/s
    margin: float = 0.5         # K below As where preheat stops
    # phase B
    Kr: float = 2.0             # 1/s   closed-loop recovery bandwidth
    Ki_r: float = 0.3           # 1/s^2
    drdt_max: float = 0.3       # 1/s   reference rate limit
    Kaw: float = 5.0            # 1/s
    T_amb_est: float = 22.0
    r_done: float = 0.99        # full-stroke latch: position cannot exceed 1, integrator would wind up
    T_limit_above_Af: float = 3.0  # K, hard over-temperature protection
    lam_scale: float = 1.0      # controller's belief about m*Lambda / true


class Controller:
    def __init__(self, p: SMAParams, cfg: CtrlCfg | None = None):
        self.p, self.c = p, cfg or CtrlCfg()
        self.map = RecoveryMap(p)
        self.T_hold = self.map.T[0] - self.c.margin
        self.Tcmd = None
        self.I_a = 0.0
        self.I_b = 0.0
        self.rref = 0.0
        self.r_est = 0.0
        self.phase = "A"
        self.done = False

    def step(self, r_target: float, T_meas: float, y_meas_mm: float) -> float:
        c, p = self.c, self.p
        if self.Tcmd is None:
            self.Tcmd = T_meas
        r_meas = float(np.clip(y_meas_mm / p.full_stroke_mm, 0.0, 1.0))
        self.r_est = max(self.r_est, r_meas)
        R = p.resistance(self.r_est)
        P_max = p.Imax ** 2 * R
        P_loss = p.h * p.surface * max(T_meas - c.T_amb_est, 0.0)
        # monotonic, rate-limited reference
        step_r = float(np.clip(r_target - self.rref, 0.0, c.drdt_max * c.Ts))
        if self.phase == "A" and T_meas >= self.T_hold - 0.1 and r_target > 0:
            self.phase = "B"; self.I_b = 0.0
        if self.phase == "A":
            goal = self.T_hold
            self.Tcmd += np.clip(c.k_soft * (goal - self.Tcmd) * c.Ts, -c.dTdt_max * c.Ts, c.dTdt_max * c.Ts)
            if goal - self.Tcmd > 0:                       # constant-rate ramp far from the target
                self.Tcmd = min(max(self.Tcmd, T_meas) , self.Tcmd + c.dTdt_max * c.Ts)
            e = self.Tcmd - T_meas
            P_un = P_loss + c.Kp * e + self.I_a
            P = float(np.clip(P_un, 0.0, P_max))
            self.I_a += c.Ts * (c.Ki * e + c.Kaw * (P - P_un))
        else:
            self.rref += step_r
            mL = p.latent_total * c.lam_scale
            e = self.rref - r_meas
            if r_target >= 1.0 and r_meas >= c.r_done:
                self.done = True
            if self.done:                                   # full stroke reached: hold temperature only
                e, self.I_b, step_r = 0.0, 0.0, 0.0
            P_un = P_loss + mL * (step_r / c.Ts + c.Kr * e + self.I_b)
            if T_meas > self.map.T[1] + c.T_limit_above_Af:
                P_un = 0.0
            P = float(np.clip(P_un, 0.0, P_max))
            self.I_b += c.Ts * (c.Ki_r * e + c.Kaw * (P - P_un) / mL)
        return float(np.sqrt(P / R))
