# MEFF recovery window: findings (2026-10-20)

## Question
Standalone-wire ANSYS (2 mm x 200 mm) recovers over ~1-1.5 K (50.8 -> 52 C). I had predicted 26 K. Model wrong or prediction wrong?

## Resolution: the earlier prediction was wrong; ANSYS is right
Source: Auricchio & Petrini 2004 (IJNME 61:807), Eq. 11 / Table I. The transformation stress carries a hardening term h*||e_tr||, and ||e_tr|| runs 0 -> eps_L, so the stress rise over the transformation is h*eps_L.

    Af - As = h * eps_L / beta = 145 * 0.04 / 5.6 = 1.04 K

Constants: h = 145 MPa, R = 167 MPa, beta = 5.6 MPa/C, T0 = 22 C, eps_L = 0.04.

    Af = T0 + R/beta                 = 51.8 C
    As = T0 + (R - h*eps_L)/beta     = 50.8 C

- ANSYS data: recovery starts ~50.8-50.95 C, complete by 52.0 C. Match.
- The earlier As = 26 C (and 35 C) is WRONG (used H/beta, missing the eps_L factor). Check thesis slide 15 and the defense deck.
- R and T0 only set WHERE the window sits, not its width.
- Ms/Mf are not defined at zero stress in this model (no strain-free martensite; Auricchio-Petrini p.813).
- Uniaxial stress shift = sqrt(3/2)*beta = 6.86 MPa/K. Loading onset = sqrt(3/2)*R = 204.5 MPa (ANSYS: 202-204.6 MPa).
- Window = (plateau stress rise)/(6.86): one h sets both plateau slope and thermal window.

## Still unexplained
- Forward transformation in ANSYS: stress rise ~2.6 MPa vs 7.1 MPa predicted; onset ~202-203 vs 204.5 MPa. Jump spans only two coarse substeps. Reverse transformation matches (1.05 K vs 1.04 K).

## Accuracy of the model vs experiment (not 100%)
- Auricchio et al. 2009 (JMEP 18:649): fitted h=100, eps_L=5.1%, beta=7.7 -> window ~0.7 K, numerical curves near-vertical; real 0.2 mm wire has steep edges plus long tails over tens of K. Model matches edges, not tails. Authors call it "good agreement" (qualitative).
- Pavia nitinol actuator springs study: Souza-Auricchio model fails below ~150 MPa; misses stroke/transformation temperatures for 1.0 mm wires; no quantitative error metrics.
- Alazzawi & Filip 2019 (Results in Materials 1:100015): one curve validated (specimen B, 128.5 g, flame heating, temperature measurement undescribed). Table 1 temperatures do not follow from its own R, H, beta (2R/beta = 60-69 K vs table Af-Ms 30-50 K). Simulated window ~5 K (exp ~8 K) vs formula 0.76 K: bending stress (0-266 MPa) shifts local temperature by up to ~39 K.
- Sources I could not open: PMC, ScienceDirect (blocked). The ANSYS text I was given is a summary, not the manual page.

## Standalone vs embedded
Stress-free As/Af are material constants; observed values shift with stress:

    As_sigma = As + sigma_start/6.86
    Af_sigma = Af + sigma_end/6.86

- Standalone free wire: 50.8 -> 51.8 C.
- In the airfoil (supervisor: As_sigma ~51, Af_sigma ~69 C): window ~18 K; implies ~118 MPa wire stress at end of recovery (derived: verify in assembly).

## Widening the window
- Alloy: raise h (cold-worked, partially annealed wire). h = W*beta/eps_L: 5 K -> ~700 MPa (plateau rise ~34 MPa), 10 K -> ~1400 MPa (~69 MPa), 20 K -> ~2800 MPa (~137 MPa, probably unrealistic). Alazzawi p.3: more dislocation -> steeper hardening slope.
- Structure (stronger lever): bias spring / airfoil stiffness raises wire stress during recovery. Extra width ~ k_b*d_eps/(sqrt(3/2)*beta).

## Controller implications
- Standalone: near-step ~1 K wide at ~51.3 C; real wire has steep edges plus tails.
- Wing: use coupled, stress-dependent window; needs airfoil stiffness at attachment as plant input.
- Tune on bench data from the real wire.

## Open items
1. Explain forward-transformation mismatch (2.6 vs 7.1 MPa).
2. Fix As in thesis text / defense deck.
3. Read end-of-recovery wire stress from the assembly.
4. Thermal/electrical data: T vs time at fixed currents, h, resistance.
5. Update plant model (corrected formulas, tunable tail term), edit requirements.md.
6. Bench data (strain-temperature at two stresses, DSC) for the actual wire.
