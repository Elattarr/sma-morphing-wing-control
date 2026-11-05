# Independent safety layer (proposal, nothing built)

The 10-check monitor shares the MCU, clock and supply with the controller, so it cannot cover a crashed or latched-up MCU.  Suggested layers, cheapest first:
1. **Hardware overcurrent latch**: INA240 output into a comparator (set about 40 A) driving a latch that forces both gate drivers off and opens the relay.  No software in the path.  The 40 A software latch in c/board is a backup only.
2. **Hardware over-temperature**: separate analog trip on a third sensor (or the MAX31856 fault pin) into the same latch.
3. **Second small MCU** (own regulator, own clock) reading T2 and the heartbeat; it holds the relay enable and the latch reset.  Heartbeat lost for > 0.5 s: relay open.
4. **Thermal fuse / 40 A fuse** in series with the wire as the last resort.
Relay opens only after the gates are low.  Part numbers and the comparator threshold are still open (hardware/pin_table.md).
