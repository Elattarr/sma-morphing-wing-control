# Bench validation plan, 0.5 mm wire (NOT DONE: no bench data exist yet)

Status: the embedded-wire model is verified (standalone wire and CFD against literature) but **not validated**.  Every controller result in this repo is model against model.  This plan is what turns it into a validation once someone runs it.

## Run
1. Setup per docs/BENCH_PROTOCOL.md (0.5 mm, 200 mm, spring 1.25 N/mm, laser range >= 10 mm, bus divider on the Vbus ADC channel).
2. Cold start, ambient measured.  Constant current 0.9 A (about 94 C, full stroke) for 145 s, then off for 250 s.  Log at >= 20 Hz: `t,I,V,T1,d` (template: data/bench_log_template.csv).  Repeat 3 times and at 2 other currents (0.6 A partial stroke, 1.1 A).
3. `python scripts/validate_bench.py log.csv --ambient <C>`

## Acceptance (fixed before the data, see script header)
RMS temperature error <= 3 K; start of stroke within 2 K; 90 % stroke within 3 K; final stroke within 10 %; V/I within 5 %.
Only the thermocouple lag is fitted.  If a criterion fails, identify h, m c, R with scripts/identify_bench.py, change the model, and re-run; do not touch the controller before the plant passes.

## Next stage, only after the plant passes
Closed-loop run with the heat-only controller at 1 target, then the 10-check monitor with injected faults (open thermocouple, stuck-on, frozen laser).

The script itself is tested only on synthetic logs (tests/test_validate_bench.py).
