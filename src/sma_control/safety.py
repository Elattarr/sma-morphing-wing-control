"""Independent safety monitor for the current-driven SMA wire.

It sits between the controller and the power stage and knows NOTHING of the controller's state: it sees the commanded current, the
MEASURED current, the control thermocouple T1, an independent over-temperature channel T2 (second thermocouple or thermal fuse with a
comparator), the position sensor d, optionally the wire voltage and a controller heartbeat. Its only output is the permitted
current and the state of a relay in series with the wire. Once tripped it latches until reset() (fail-safe state: power OFF; the bias
spring then pulls the wire to the extended end, which is the mechanically safe end).

Checks (each one targets a failure the controller cannot see by itself):
  1. HARD OVER-TEMPERATURE   max(T1, T2) > T_trip                                      -> 'over_temperature'
  2. CHANNEL DISAGREEMENT    |T1 - T2| > T_disagree (beyond T_disagree_dyn while T changes fast: sensor lags)  -> 'tc_disagree'
  3. THERMOCOUPLE ENVELOPE   energy-balance bounds from the measured power; T1 must lie between them
     (heat capacity m c between 0.3x and 2x nominal, latent heat up to 3x nominal, h between h_min and h_max;
     two observers, small and large m c, because heating and cooling need opposite extremes):
         upper:  m c dT = I^2 R_max - h_min S (T - T_amb_hi)                   (no latent heat credit)
         lower:  m c dT = I^2 R_min - h_max S (T - T_amb_lo) - latent sink      (full convection, latent heat absorbed first)
     catches an open or shorted thermocouple, a sensor stuck at a plausible value while the heater runs, a wrong channel.
  4. STALL                   T1 above A_f,max + stall_dT for t_stall of heating (I_meas > 0.3 A) without any motion of the position sensor
         (jammed wire, stuck/broken position sensor): needs the stroke from homing                -> 'stall'
  5. CURRENT MISMATCH        measured current > commanded + tolerance (power stage stuck on)                   -> 'current_mismatch'
  6. WIRE RESISTANCE         V/I outside [R_min, R_max] (broken wire, loose contact)                           -> 'resistance'
  7. WATCHDOG                heartbeat counter not changing (controller hung, current command frozen)          -> 'watchdog'
  9. FROZEN THERMOCOUPLE     T1 bit-identical for t_tc_freeze while heated (I_meas > 0.3 A) or hot (> T_amb + 5 K): a real reading always has noise    -> 'tc_frozen'
  8. SENSOR RANGE            T1 or T2 outside [-20, 250] C (open circuit rails)                                -> 'sensor_range'

[ASSUMED] bounds: convection coefficient h in [h_min, h_max] = [5, 100] W/m2/K, ambient within +/- 10 K of the value at arming.
Nothing here is validated on hardware; the thresholds are chosen so that no check fires in the evaluation sequences of
scripts/evaluate_controller.py (zero false trips) and each injected fault is caught before the wire passes T_trip.
"""
from dataclasses import dataclass
import numpy as np
from .params import SMAParams


@dataclass
class SafetyCfg:
    T_trip: float = 130.0            # C   hard trip, above the controller's absolute limit (125 C) so it never acts in normal use
    T_disagree: float = 12.0         # K   static disagreement between T1 and T2
    T_disagree_dyn: float = 25.0     # K   allowed while |dT/dt| of either channel exceeds 3 K/s (thermocouple lags differ)
    h_min: float = 5.0               # W/m2/K  [ASSUMED] lowest credible convection (still air, enclosure)
    h_max: float = 100.0             # W/m2/K  [ASSUMED] highest credible convection (forced air)
    dT_amb: float = 10.0             # K   ambient uncertainty after arming
    band: float = 6.0                # K   envelope margin (noise, latent heat released on cooling)
    T_lat_lo: float = 55.0           # C   below this the lower bound takes no latent-heat sink (lowest credible A_s)
    stall_dT: float = 5.0            # K   above A_f,max
    t_stall: float = 8.0             # s
    stall_frac: float = 0.15         # of the homing stroke: smaller motion counts as none
    I_tol_abs: float = 0.4           # A   measured-current tolerance
    I_tol_rel: float = 0.5
    t_I_mismatch: float = 1.0        # s
    R_tol: float = 0.45              # relative tolerance on the wire resistance
    t_tc_freeze: float = 10.0        # s   unchanged T1 reading while the wire is heated or hot
    n_freeze: int = 10               # samples with an unchanged heartbeat
    mc_min_frac: float = 0.3         # [ASSUMED] smallest credible m*c as a fraction of the nominal value (c as low as 250 J/kg/K): upper bound
    mc_max_frac: float = 2.0         # [ASSUMED] largest credible m*c (clamps, wires, glue): lower bound
    lat_max_frac: float = 3.0        # [ASSUMED] largest credible latent heat as a multiple of the nominal one: lower bound
    tc_lag_max: float = 4.0          # s   [ASSUMED] slowest credible thermocouple time constant; the lower bound is lagged by it
    T_range: tuple = (-20.0, 250.0)


class SafetyMonitor:
    def __init__(self, p: SMAParams, Ts: float, cfg: SafetyCfg | None = None):
        self.p, self.Ts, self.c = p, Ts, cfg or SafetyCfg()
        self.mc, self.mL, self.S = p.heat_capacity, p.latent_total, p.surface
        self.mcs, self.mL_max = (self.mc * self.c.mc_min_frac, self.mc * self.c.mc_max_frac), self.mL * self.c.lat_max_frac
        Rs = [p.resistance(r) for r in np.linspace(0.0, 1.0, 11)]
        self.R_min, self.R_max = min(Rs), max(Rs)
        self.reset()

    # ------------------------------------------------------------------ state
    def reset(self):
        self.tripped, self.reason, self.t = False, None, 0.0
        self.armed = False
        self.T_lo = self.T_lo_f = self.T_hi = self.Tamb0 = None
        self.To = None
        self._I_prev, self._t_Imis, self._t_hot, self._d_hot0 = 0.0, 0.0, 0.0, None
        self._beat, self._n_same, self._T1_prev, self._T2_prev = None, 0, None, None
        self._T1_frz, self._t_frz = None, 0.0
        self.stroke_mm, self.Af_max = None, None

    def arm(self, T1: float, T2: float | None = None):
        """call once at power-up with the wire cold and the heater off"""
        lo, hi = self.c.T_range
        if not (lo < T1 < hi) or (T2 is not None and not (lo < T2 < hi)):
            self._trip("sensor_range")
            return
        if T2 is not None and abs(T1 - T2) > self.c.T_disagree:
            self._trip("tc_disagree")
            return
        self.Tamb0 = T1
        self.T_lo, self.T_hi = T1 - self.c.dT_amb, T1 + self.c.dT_amb
        # two observers (small and large m c): heating needs the small one for the upper bound and the large one for the lower bound,
        # cooling the other way round; the bounds are the max / min over both
        self.To = [dict(hi=self.T_hi, lo=self.T_lo, lof=self.T_lo, lat=self.mL_max) for _ in self.mcs]
        self.T_lo_f = self.T_lo
        self.armed = True

    def set_calibration(self, stroke_mm: float, Af_max: float):
        """after homing: stroke [mm] between the two end points and the highest end-of-recovery temperature"""
        self.stroke_mm, self.Af_max = stroke_mm, Af_max

    def _trip(self, why: str):
        if not self.tripped:
            self.tripped, self.reason = True, why

    # ------------------------------------------------------------------ one sample
    def step(self, I_cmd: float, I_meas: float, T1: float, T2: float | None = None, d: float | None = None,
             V: float | None = None, beat: int | None = None):
        """returns (I_permitted, relay_closed). Call every Ts with the newest samples."""
        c, dt = self.c, self.Ts
        self.t += dt
        if not self.armed and not self.tripped:
            self.arm(T1, T2)
        if self.tripped:
            return 0.0, False
        Tm = max(T1, T2) if T2 is not None else T1
        lo, hi = c.T_range

        # 8 + 1: range and hard over-temperature
        if not (lo < T1 < hi) or (T2 is not None and not (lo < T2 < hi)):
            self._trip("sensor_range")
        if Tm > c.T_trip:
            self._trip("over_temperature")

        # 2: channel disagreement (allowed larger while temperature changes fast)
        if T2 is not None:
            fast = False
            if self._T1_prev is not None:
                fast = max(abs(T1 - self._T1_prev), abs(T2 - self._T2_prev)) / dt > 3.0
            if abs(T1 - T2) > (c.T_disagree_dyn if fast else c.T_disagree):
                self._trip("tc_disagree")
        self._T1_prev, self._T2_prev = T1, T2 if T2 is not None else T1

        # 3: thermocouple envelope from the MEASURED power
        P_hi, P_lo = I_meas * I_meas * self.R_max, I_meas * I_meas * self.R_min
        Tamb_hi, Tamb_lo = self.Tamb0 + c.dT_amb, self.Tamb0 - c.dT_amb
        for o, mc in zip(self.To, self.mcs):
            o["hi"] += dt * (P_hi - c.h_min * self.S * (o["hi"] - Tamb_hi)) / mc
            T_old = o["lo"]
            T_new = T_old + dt * (P_lo - c.h_max * self.S * (T_old - Tamb_lo)) / mc
            if T_new > c.T_lat_lo and T_new > T_old and o["lat"] > 0.0:           # heating through the transformation: latent sink
                dE = mc * (T_new - max(T_old, c.T_lat_lo))
                take = min(dE, o["lat"])
                o["lat"] -= take
                T_new -= take / mc
            elif T_new < c.T_lat_lo - 5.0:
                o["lat"] = self.mL_max
            o["lo"] = T_new
            o["lof"] += dt / c.tc_lag_max * (o["lo"] - o["lof"])                   # a lagging sensor trails the true temperature
        self.T_hi = max(o["hi"] for o in self.To)
        self.T_lo = min(o["lo"] for o in self.To)
        self.T_lo_f = min(o["lof"] for o in self.To)
        if T1 > self.T_hi + c.band:
            self._trip("tc_reads_high")
        if T1 < self.T_lo_f - c.band:
            self._trip("tc_reads_low")

        # 4: stall (needs the homing stroke)
        if self.stroke_mm is not None and d is not None:
            if T1 > self.Af_max + c.stall_dT:
                if self._d_hot0 is None:
                    self._d_hot0, self._t_hot = d, 0.0
                if I_meas > 0.3:                                    # only a HEATED wire that does not move is a stall (a hot wire cooling at rest is not)
                    self._t_hot += dt
                if abs(d - self._d_hot0) > c.stall_frac * self.stroke_mm:
                    self._d_hot0, self._t_hot = d, 0.0
                if self._t_hot > c.t_stall:
                    self._trip("stall")
            else:
                self._d_hot0, self._t_hot = None, 0.0

        # 9: frozen thermocouple (exactly the same reading for t_tc_freeze while the wire is heated or hot)
        if self._T1_frz is not None and T1 == self._T1_frz and (I_meas > 0.3 or T1 > self.Tamb0 + 5.0):
            self._t_frz += dt
            if self._t_frz > c.t_tc_freeze:
                self._trip("tc_frozen")
        else:
            self._t_frz = 0.0
        self._T1_frz = T1

        # 5: current mismatch
        if I_meas > self._I_prev * (1.0 + c.I_tol_rel) + c.I_tol_abs:
            self._t_Imis += dt
            if self._t_Imis > c.t_I_mismatch:
                self._trip("current_mismatch")
        else:
            self._t_Imis = 0.0

        # 6: wire resistance
        if V is not None and I_meas > 0.3:
            R = V / I_meas
            if not (self.R_min * (1.0 - c.R_tol) < R < self.R_max * (1.0 + c.R_tol)):
                self._trip("resistance")

        # 7: watchdog
        if beat is not None:
            self._n_same = self._n_same + 1 if beat == self._beat else 0
            self._beat = beat
            if self._n_same >= c.n_freeze:
                self._trip("watchdog")

        if self.tripped:
            return 0.0, False
        self._I_prev = float(min(max(I_cmd, 0.0), self.p.Imax))
        return self._I_prev, True
