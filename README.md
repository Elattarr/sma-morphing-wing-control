# Standalone SMA wire control

Controller and plant simulation for a free one-way SMA wire (MEFF model, Alazzawi & Filip 2019, specimen C) pulled back by a bias spring and driven by current.

Simulation-verified only; nothing is validated on a real wire. See `docs/STATUS.md` (done / not done / limits), `docs/ASSUMPTIONS.md` (every assumed number) and `docs/BENCH_PROTOCOL.md` (hardware spec and commissioning gates).

- `src/sma_control/`: plant, bias spring, parallel-element plant, cycle controller with homing, safety monitor, identification
- `c/`: portable C99 port (controller, homing, monitor), replay-tested against Python
- `scripts/`: evaluations (`evaluate_controller.py`, `evaluate_materials.py`, `identify_bench.py`)
- `docs/eval_*.txt`: evaluation results; `docs/figs_*`: figures

Tests: `python -m pytest`
