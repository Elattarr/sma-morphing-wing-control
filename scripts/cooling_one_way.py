"""ONE-WAY wire (pure MEFF): cooling at zero stress does nothing; a bias stress brings the transformation back.
python scripts/cooling_one_way.py [figure_dir]"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src")); sys.path.insert(0, str(Path(__file__).parent))
from sma_control import SMAParams
from sma_control.cycle_plant import SMAWireCycle
from cooling_common import cycle, first_forward, style, time_panel, T_AMB

SIGMAS = (0.0, 100.0, 204.0, 250.0)


def run():
    out = []
    for s in SIGMAS:
        w = SMAWireCycle(SMAParams(), T0=T_AMB, xi0=1.0)
        L, nh = cycle(w, s)
        out.append(dict(sigma=s, w=w, L=L, nh=nh, A_f=w.u.A_f(s), T_fwd=first_forward(L, nh), T_theory=w.T_forward_start(s),
                        d_min=L[:, 3].min(), d_end=L[-1, 3], xi_end=L[-1, 2], T_max=L[:, 1].max()))
    return out


def report(out):
    print("ONE-WAY wire, pure MEFF. Start xi = 1 at 22 C, heat 15 A to A_f+3 K, cool 900 s with I = 0.")
    print(f"{'sigma [MPa]':>11s} {'A_f [C]':>8s} {'Tmax':>6s} {'d hot [mm]':>11s} {'fwd start sim':>14s} {'theory':>7s} {'d end [mm]':>11s} {'xi end':>7s}")
    for o in out:
        print(f"{o['sigma']:11.0f} {o['A_f']:8.1f} {o['T_max']:6.1f} {o['d_min']:11.3f} {o['T_fwd']:14.2f} {o['T_theory']:7.2f} {o['d_end']:11.3f} {o['xi_end']:7.3f}")
    print("fwd start sim = nan: no re-extension down to 22 C.  d includes the elastic strain sigma/E under bias.")


def figure(out, outdir):
    plt = style()
    cols = ["#52514e", "#2a78d6", "#1baf7a", "#eb6834"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for o, c in zip(out, cols):
        ax[0].plot(o["L"][:, 1], o["L"][:, 3], color=c, lw=1.8, label=f"bias {o['sigma']:.0f} MPa")
    ax[0].set_xlabel("wire temperature [°C]"); ax[0].set_ylabel("displacement [mm]"); ax[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    ax[0].set_title("One-way: displacement against temperature", loc="left", fontsize=10)
    time_panel(ax[1], out[3]["L"], "250 MPa bias against time: black = displacement, orange = T")
    fig.tight_layout(); Path(outdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(outdir) / "fig8a_one_way_cooling.png", dpi=160, facecolor="#fcfcfb")


if __name__ == "__main__":
    o = run(); report(o)
    if len(sys.argv) > 1:
        figure(o, sys.argv[1])
