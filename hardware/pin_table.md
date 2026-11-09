# SMA drive: pin-by-pin connection table (design, not built)

"Checked" = read in a datasheet this session. "Not checked" = my pin/value from memory or a placeholder class, verify before ordering.

## Power path (35 A)
| From | To | Note | Status |
|---|---|---|---|
| Supply + (5 V, 40 A) | F1 fuse 40 A | fast blade, cable ≥10 mm² (35 A) | not checked |
| F1 | K1 contact | contact ≥40 A DC | not checked |
| K1 | +V rail: C1, Q2 drain, L1 | C1 2×1000 µF polymer + ceramics at the FETs | not checked |
| L1 (~10 µH, ≥60 A sat) | NiTi wire | wire: model 0.051-0.064 Ω; contact R unmeasured (the earlier "≈0.1 Ω" was a guess) | not checked |
| Wire | SW node | | |
| SW node | Q1 drain and Q2 source | Q2 = freewheel, drain on +V | |
| Q1 source | R1 shunt 1 mΩ → GND star | 4-terminal, ≥5 W | not checked |
| Supply − | GND star | one point only | |

## MOSFETs Q1, Q2: BSC014N04LS (PG-TDSON-8), checked
VDS 40 V, RDS(on) 1.1 mΩ typ / 1.4 max at 10 V, QG 61 nC typ, VGS(th) 1.2–2 V, Rth(j-a) 50 K/W on 6 cm² Cu, body diode 0.82 V typ. Pinout not read: take it from the package drawing.

## Gate driver U1: UCC27211D (SOIC-8), checked
VDD 8–17 V, inputs TTL (VIH 2.3 V, so 3.3 V logic works; the UCC27210 is the CMOS one and would not), peak 3.7 A source / 4.5 A sink min, HB 120 V max, UVLO 7.0 V rising. No built-in dead time, so HRTIM must make it.
| Pin | Name | Connects to |
|---|---|---|
| 1 | VDD | +12 V, C3 1 µF to VSS, D1 anode |
| 2 | HB | D1 cathode, C2 1 µF to HS |
| 3 | HO | Rg 2.2 Ω → Q2 gate |
| 4 | HS | SW node, Q2 source |
| 5 | HI | PA9 (HRTIM1_CHA2) |
| 6 | LI | PA8 (HRTIM1_CHA1) |
| 7 | VSS | Q1 source (Kelvin), 12 V return |
| 8 | LO | Rg 2.2 Ω → Q1 gate |

## Current amplifier U2: INA240A2 (SOIC-8)
Checked: gain 50 V/V (A2), common mode −4…80 V, offset ±25 µV max, gain error 0.2 % max. Not checked: pin numbers, supply range, output swing to GND. IN+ = shunt FET side, IN− = shunt GND side, REF1 and REF2 to GND (if swing to GND is poor near 0 A, bias REF), OUT → 100 Ω / 1 nF → PA0. 35 A → 35 mV → 1.75 V.

## Thermocouple amplifiers U3, U4: MAX31856 (TSSOP-14)
Checked: types K J N E T R S B, 19-bit, ±0.15 % full scale, cold junction ±0.7 °C, inputs ±45 V. Not checked: pin numbers, SPI mode, supply (use a breakout board first). SCK PB13, MISO PB14, MOSI PB15 shared; CS U3 = PB12, U4 = PB1. Junctions electrically insulated from the wire.

## MCU: NUCLEO-G474RE (pins not checked against the board manual)
| Pin | Function | To |
|---|---|---|
| PA8 | HRTIM1_CHA1 (AF13) | U1 pin 6 |
| PA9 | HRTIM1_CHA2 (AF13) | U1 pin 5 |
| PA0 | ADC1_IN1 | U2 OUT |
| PA1 | ADC1_IN2 | laser 4–20 mA across 150 Ω, RC 1 kΩ / 10 nF |
| PB13/14/15 | SPI2 | U3, U4 |
| PB12, PB1 | GPIO | CS |
| PB0 | GPIO | Q3 2N7002 gate → K1 coil (12 V, flyback 1N4007) |
| PA2/PA3 | LPUART1 | ST-LINK USB to laptop |

## Rules that are part of the design
1. Relay closes only with LO low and HI low; opens first on any fault.
2. HRTIM dead time ≥ 100 ns (placeholder, set from the QG and gate-resistor numbers).
3. ADC triggered by HRTIM at the middle of the Q1 on-time.
4. Thermocouple junctions never touch the wire or its clamps.
5. Laser sees the wire end (0–6.64 mm), not the trailing edge.

## Still open
Second safety processor, return mechanism, real shunt and relay part numbers, SPICE run of the power stage.

## Datasheet check 2026-11-05
Checked: PA0 ADC1_IN1, SPI2 PB12-15, LPUART1 PA2/PA3, INA240 VS 2.7-5.5 V and output swing, MAX31856 registers.
Unchecked: INA240 / MAX31856 pin numbers, HRTIM pins, PA1, Vbus divider pin, SPI mode, relay, fuse, shunt, laser (range >= 10 mm needed).

Laser: range 10-20 mm required (stroke 6.6 mm); 12-bit ADC over 4-20 mA gives about 0.005 mm per count at 15 mm, a 100 mm range gives 0.034 mm and fails the noise requirement.
