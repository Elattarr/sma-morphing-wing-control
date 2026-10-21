"""Step 2: heating model of the free SMA wire (current -> wire temperature).

Equation (lumped energy balance, same form as the thesis):
    rho*C*V*dT/dt + lam*V*dxi/dt + h*S*(T - Tinf) = I^2 * R(xi)
      left  : stored heat + latent heat (reverse transformation) + heat lost to air
      right : Joule heating
    xi  = austenite fraction (0 = martensite, 1 = austenite), rises with T between As and Af
    R   = resistance, mix of martensite and austenite values
PLACEHOLDER: xi(T) is a straight ramp between As=26 C and Af=52 C until we have the ANSYS curve.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "data"))
from thesis_constants import WIRE, T_AMBIENT, AS_FREE, AF_FREE, FROZEN_STRAIN

d, L = WIRE["d_mm"] * 1e-3, WIRE["L_mm"] * 1e-3
A = np.pi * d**2 / 4          # cross-section [m2]
V = A * L                     # volume [m3]
S = np.pi * d * L             # lateral surface [m2]
rho, C = WIRE["rho"], WIRE["C"]
lam = WIRE["latent_J_per_mm3"] * 1e9      # J/m3
h = WIRE["h"]

def xi_of_T(T, T1=AS_FREE, T2=AF_FREE):
    """austenite fraction vs temperature (placeholder ramp)"""
    return float(np.clip((T - T1) / (T2 - T1), 0.0, 1.0))

def resistance(xi):
    return ((xi * WIRE["rho_E_A"] + (1 - xi) * WIRE["rho_E_M"]) * L) / A

def heat(I, t_end, T_inf=T_AMBIENT, h_=h, T1=AS_FREE, T2=AF_FREE, T_stop=None):
    """simulate constant current I [A]. returns t, T, xi, energy [J]"""
    def f(t, y):
        T = y[0]
        xi = xi_of_T(T, T1, T2)
        dxi_dT = 1.0 / (T2 - T1) if T1 < T < T2 else 0.0
        P = I**2 * resistance(xi)
        return [(P - h_ * S * (T - T_inf)) / (rho * C * V + lam * V * dxi_dT)]
    events = []
    if T_stop is not None:
        ev = lambda t, y: y[0] - T_stop
        ev.terminal = True
        events = [ev]
    s = solve_ivp(f, [0, t_end], [T_inf], events=events, max_step=0.01)
    xi = np.array([xi_of_T(T, T1, T2) for T in s.y[0]])
    E = np.trapezoid(I**2 * resistance(xi), s.t)
    return s.t, s.y[0], xi, E

def current_for_time(t_target, T_final, **kw):
    """constant current that reaches T_final in t_target seconds"""
    def miss(I):
        t, T, xi, E = heat(I, 4 * t_target + 20, T_stop=T_final, **kw)
        return (t[-1] if T[-1] >= T_final - 1e-6 else 1e3) - t_target
    return brentq(miss, 3.0, 80.0, xtol=1e-3)

if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    T_final = AF_FREE   # heat the free wire to Af
    print(f"Free wire: 22 C -> {T_final:.0f} C")
    print("sensible heat [J]:", round(rho*C*V*(T_final-T_AMBIENT), 1), "| latent [J]:", round(lam*V, 1))
    print("hold power at Af [W]:", round(h*S*(T_final-T_AMBIENT), 2))
    fig, ax = plt.subplots(figsize=(6, 4))
    for tt in [2.5, 5.0, 10.0]:
        I = current_for_time(tt, T_final)
        t, T, xi, E = heat(I, tt + 5, T_stop=T_final)
        print(f"t={tt:>4}s  I={I:5.1f} A  J={I/(A*1e6):4.1f} A/mm2  E={E:6.1f} J  meanP={E/t[-1]:5.1f} W")
        ax.plot(t, T, label=f"{I:.1f} A")
    ax.set_xlabel("time [s]"); ax.set_ylabel("wire temperature [C]"); ax.grid(alpha=.3); ax.legend(title="current")
    plt.tight_layout(); plt.savefig(root / "figs" / "heating_free_wire.png", dpi=140)
