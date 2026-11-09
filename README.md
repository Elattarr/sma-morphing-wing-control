# Standalone SMA wire control

Controller and plant simulation for NiTi shape-memory-alloy (SMA) wire actuation, written for a morphing-wing master thesis (NACA 2412, 1000 mm chord, one 2 mm × 200 mm wire embedded in the airfoil).

**Status: simulation-verified only. Nothing is validated on a real wire or wing.** Everything in `docs/results/eval_*.txt` is model against model, or model against ANSYS (itself a simulation). Read `docs/STATUS.md` for what is done, what failed and the limits before you trust any number.

## Start here

| If you want to... | Open |
|---|---|
| know what is done, what failed, what is unproven | `docs/STATUS.md` |
| see the drive circuit working (animated, with fault injection) | `hardware/drive_circuit_animation.html` (open in a browser) |
| understand the safety design | `docs/SAFETY_PROCESSOR.md`, `src/sma_control/safety.py` |
| see the controller results animated and charted | `docs/heatonly_animation.html` (browser), `docs/figures/heatonly_robust_bars.png`, `docs/figures/heatonly_windows.png` |
| validate the model on a real 0.5 mm wire | `docs/BENCH_PROTOCOL.md`, `docs/BENCH_VALIDATION_PLAN.md`, `scripts/validate_bench.py` |

## Repository layout

```
src/sma_control/   Python: wire model, controller, safety monitor, trailing-edge map, hardware interface
c/                 C99 port of controller + monitor (replay-tested against Python)
c/board/           board layer: 20 kHz current loop, slow loop, relay logic, simulated power stage (PC only)
tests/             pytest suite (about 150 tests, about 2 minutes)
scripts/           evaluations, C replay export, board chain run, bench identification and validation
data/              ANSYS exports used for fitting, bench log template
docs/              STATUS (read first), ASSUMPTIONS, REQUIREMENTS, bench docs, safety processor, MEFF findings
docs/results/      printed outputs of the evaluations
docs/figures/      result plots
hardware/          circuit animation, wiring diagrams, pin table (design only, nothing built)
```

## Install and test

```bash
git clone https://github.com/Elattarr/sma-morphing-wing-control.git
cd sma-morphing-wing-control
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m pytest                                        # about 2 minutes; C-port tests need gcc and are skipped without it
```

Python 3.10 or newer (the code uses `X | None` annotations).

## Run the heat-only evaluation

The thesis wire is selected by environment variables inside `scripts/evaluate_heatonly.py` (2 mm wire, 35 A, 6.36 N/mm structure stiffness). Run from the repo root:

```bash
PYTHONPATH=src python scripts/evaluate_heatonly.py steps     # step responses: delay, rise, settling, overshoot
PYTHONPATH=src python scripts/evaluate_heatonly.py ramps
PYTHONPATH=src python scripts/evaluate_heatonly.py disturb
PYTHONPATH=src python scripts/evaluate_heatonly.py robust    # 16 parameter-error cases
PYTHONPATH=src python scripts/evaluate_heatonly.py sensors
PYTHONPATH=src python scripts/evaluate_heatonly.py safety    # 30 injected faults
```

Each section prints a table; the committed outputs are `docs/results/eval_heatonly_*.txt`.

## Use the heat-only controller

```python
from sma_control.params import SMAParams
from sma_control.cycle_controller import CycleController, cfg_heat_only
from sma_control.spring import BiasSpring

P = SMAParams(d=2e-3, Imax=35.0)                     # 2 mm NiTi wire, 35 A limit
K_structure = 6.36                                    # N/mm at the wire (from ANSYS), not a bias spring
ctl = CycleController(P, 2.2, cfg_heat_only(P, Ts=0.05), K_spring=BiasSpring(2.2, K_structure).K(P))
# every sample: I = ctl.step(r_target, T_wire_degC, wire_end_displacement_mm)
```

`hardware_io.run_step` wraps this with the safety monitor and the TE map. See `tests/test_heat_only_unit.py` for a complete simulated loop.

## Results in short (simulation)

- Wide windows (18 K, 30 K, 25 K with A_f 73 °C) reach 4 of 4 targets in simulation, but the current chatters (0 to about 30 A) in the last strokes: not tuned yet.
- Heat-only on the 18 K wing window: 98.6 % of the stroke at 68.8 °C peak (not 100 %).
- Full stroke: settles in 14.5 s, overshoot below 1 %, steady-state error below 0.4 %, wire peak 54–57 °C, about 220 J.
- 16 of 18 parameter-error cases reach all four targets (0.25 / 0.5 / 0.75 / 1.0 of the stroke), including a plant fitted to the ANSYS recovery curve.
- Needs displacement noise of about 0.03 mm or better on the wire stroke; fails at 0.1 mm.
- All 30 injected safety faults end safe.

## Limits (not hidden)

- No hardware validation; gains and safety thresholds are starting points for a bench.
- Heating only. The wing return path is not controlled or modelled.
- Thermal model is air convection, not wire embedded in CFRP.
- The TE map and recovery window come from ANSYS simulations.
- Two stress tests fail (5 K window spread; 25 % tail over 20 K); both are much wider than the ANSYS data.

## Reference

Alazzawi & Filip (2019), MEFF model for NiTi, specimen C. Model constants and assumptions: `docs/ASSUMPTIONS.md`.

## License

MIT, see `LICENSE`. Use at your own risk: this code can drive a heater at 35 A; do not run it on real hardware without the commissioning gates in `docs/BENCH_PROTOCOL.md` and an independent hardware over-temperature cutout.
