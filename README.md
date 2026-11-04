# Standalone SMA wire control

Controller and plant simulation for NiTi shape-memory-alloy (SMA) wire actuation, written for a morphing-wing master thesis (NACA 2412, 1000 mm chord, one 2 mm × 200 mm wire embedded in the airfoil).

**Status: simulation-verified only. Nothing is validated on a real wire or wing.** Everything in `docs/results/eval_*.txt` is model against model, or model against ANSYS (itself a simulation). Read `docs/STATUS.md` for what is done, what failed and the limits before you trust any number.

## Repository layout

```
src/sma_control/   the library: plant, controller, safety monitor, TE map, hardware interface
c/                 C99 port of the controller and monitor, replay-tested against Python
tests/             114 tests
scripts/           evaluations (evaluate_heatonly.py is the main one), C-replay export, bench identification
data/              ANSYS exports used for fitting
docs/              STATUS.md (read first), ASSUMPTIONS, BENCH_PROTOCOL, REQUIREMENTS, MEFF_window_findings
docs/results/      printed outputs of the evaluations
docs/figures/      the heat-only result plots and ANSYS summaries
hardware/          wiring diagram, pin table, drive circuit (design only, nothing built)
```

| Part | What it is |
|---|---|
| `src/sma_control/cycle_controller.py` | Two controllers in one class: the bias-spring cycle controller with homing, and the **heat-only (one-way) controller** for the embedded wire (`cfg_heat_only`): temperature creep into the transformation window, then a slow position loop; no spring, no return; first-stroke learner; fault states |
| `src/sma_control/safety.py` | Independent monitor between controller and power stage, 10 checks: over-temperature, thermocouple disagreement, energy-balance envelope, stall, current mismatch, resistance, watchdog, sensor range, frozen thermocouple, frozen position sensor |
| `src/sma_control/te_map.py` | Wire recovery r to trailing-edge deflection in mm, from the ANSYS Up 40 mm run (concave, not linear) |
| `src/sma_control/hardware_io.py` | `RigIO` interface and one control-loop step; you implement ADC / PWM / relay in a subclass (untested on hardware) |

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
