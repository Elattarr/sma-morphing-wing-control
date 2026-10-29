"""TWO-WAY wire (assumed extension, NOT MEFF): cooling at zero stress returns part of the stored displacement.
python scripts/cooling_two_way.py [figure_dir]"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src")); sys.path.insert(0, str(Path(__file__).parent))
from sma_control import SMAParams
from sma_control.two_way import SMAWireTwoWay, TwoWayParams
from cooling_common import cycle, first_forward, style, time_panel, T_AMB

CASES = [TwoWayParams(35.0, 0.3), TwoWayParams(30.0, 0.3), TwoWayParams(45.0, 0.3), TwoWayParams(35.0, 0.6)]


def run():
    out = []
    for tw in CASES:
        w = SMAWireTwoWay(SMAParams(), T0=T_AMB, xi0=1.0, tw=tw)
        L, nh = cycle(w, 0.0)
        out.append(dict(tw=tw, w=w, L=L, nh=nh, T_fwd=first_forward(L, nh), d_end=L[-1, 3], xi_end=L[-1, 2],
                        stored=w.p.L * 1e3 * w.u.eps_tr))
    return out


def report(out):
    print("TWO-WAY wire (assumed extension). sigma = 0. Start xi = 1 at 22 C, heat 15 A to A_f+3 K, cool 900 s with I = 0.")
    print(f"{'Ms [C]':>7s} {'xi_cap':>7s} {'sigma_int [MPa]':>16s} {'fwd start sim':>14s} {'d end [mm]':>11s} {'xi end':>7s} {'stored [mm]':>12s}")
    for o in out:
        tw = o["tw"]
        print(f"{tw.Ms:7.1f} {tw.xi_cap:7.2f} {o['w'].sig_int:16.1f} {o['T_fwd']:14.2f} {o['d_end']:11.3f} {o['xi_end']:7.3f} {o['stored']:12.3f}")
    print("sigma_int = stress that puts the zero-stress forward start at Ms (MEFF numbers: large because the zero-stress")
    print("hysteresis of this parameter set is about 60 K). Ms and xi_cap are assumptions.")


def figure(out, outdir):
    plt = style()
    cols = ["#7a4fd1", "#2a78d6", "#1baf7a", "#eb6834"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for o, c in zip(out, cols):
        ax[0].plot(o["L"][:, 1], o["L"][:, 3], color=c, lw=1.8, label=f"Ms {o['tw'].Ms:.0f} °C, cap {o['tw'].xi_cap:.1f}")
    ax[0].set_xlabel("wire temperature [°C]"); ax[0].set_ylabel("displacement [mm]"); ax[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    ax[0].set_title("Two-way (assumed), σ = 0: displacement against temperature", loc="left", fontsize=10)
    time_panel(ax[1], out[0]["L"], "Ms 35 °C, cap 0.3 against time: black = displacement, orange = T")
    fig.tight_layout(); Path(outdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(outdir) / "fig8b_two_way_cooling.png", dpi=160, facecolor="#fcfcfb")


if __name__ == "__main__":
    o = run(); report(o)
    if len(sys.argv) > 1:
        figure(o, sys.argv[1])
