"""Hardware interface for the control loop: what the real rig has to provide, and a simulated implementation.

The control loop (`run_step`) only talks to `RigIO`; the bench code subclasses it (ADC / DAC / relay). Nothing here is tested on hardware.
Units: temperatures in C, displacement of the wire in mm (positive as in the model: extended wire = start value), current in A.
Order per sample: read sensors -> controller -> safety monitor -> write the PERMITTED current and the relay. The monitor sees the
MEASURED current, the second thermocouple and a heartbeat counter that the loop increments every sample.
"""
from dataclasses import dataclass
from typing import Optional
from .cycle_controller import CycleController
from .safety import SafetyMonitor
from .te_map import TEMap


class RigIO:
    """override in the bench code"""
    def read_T1(self) -> float: raise NotImplementedError          # control thermocouple
    def read_T2(self) -> Optional[float]: return None               # independent over-temperature channel
    def read_disp_mm(self) -> float: raise NotImplementedError      # displacement of the wire (end block), NOT the trailing edge
    def read_I_meas(self) -> float: raise NotImplementedError       # measured wire current
    def read_V(self) -> Optional[float]: return None                # wire voltage (resistance check)
    def write_current(self, amps: float) -> None: raise NotImplementedError
    def write_relay(self, closed: bool) -> None: raise NotImplementedError


@dataclass
class LoopState:
    beat: int = 0
    I_cmd: float = 0.0
    I_perm: float = 0.0
    relay: bool = False
    fault: Optional[str] = None          # controller fault or monitor trip reason


def run_step(io: RigIO, ctl: CycleController, mon: SafetyMonitor, te_target_mm: float, te_map: TEMap, st: LoopState) -> LoopState:
    """one sample of the loop; returns the updated LoopState"""
    T1, T2, d = io.read_T1(), io.read_T2(), io.read_disp_mm()
    I_meas, V = io.read_I_meas(), io.read_V()
    r_target = te_map.te_to_r(te_target_mm)
    st.I_cmd = ctl.step(r_target, T1, d)
    st.beat += 1
    st.I_perm, st.relay = mon.step(st.I_cmd, I_meas, T1, T2, d, V, st.beat)
    io.write_current(st.I_perm)
    io.write_relay(st.relay)
    st.fault = mon.reason if mon.tripped else ctl.fault
    return st
