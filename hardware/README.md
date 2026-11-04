# Hardware concept (design only, nothing built)

5 V / 35 A drive for the 2 mm × 200 mm NiTi wire: half-bridge with a freewheel MOSFET, ground-side shunt, STM32G474 timer PWM. Sized by calculation and simulation; no SPICE run, no bench result.

| File | What it is |
|---|---|
| `wiring_diagram.png` / `.svg` | Schematic-style wiring with real part numbers |
| `real_components.png` | The same connections drawn as physical parts |
| `pin_table.md` | Pin-by-pin table; says which values were read in datasheets and which are placeholders |
| `drive_circuit_animation.html` | Open in a browser: animated circuit with fault injection (model data) |

Checked in datasheets: UCC27211 (pins, TTL thresholds, 8–17 V supply), BSC014N04LS, INA240A2 gain and range, MAX31856 type. Not checked: fuse, relay, shunt, laser, STM32 pin numbers. Open: wire-voltage input for the resistance check, second safety processor, return mechanism, board layer and fast current loop.
