# Plan
1. Get standalone ANSYS data (temperature, strain) -> data/
2. Output map: temperature -> contraction (fit or interpolate), error vs ANSYS
3. Heating model (energy balance), check against thesis energies 173/182/199 J
4. PI controller, then PI + feed-forward
5. Tests (pytest) for R1-R7, then README
