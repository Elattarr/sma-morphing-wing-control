"""MEFF parameter study (see docs: MEFF_Terms_Referenced_Study). Needs numpy + scipy. Run: python scripts/meff_study.py
Cited relations are tagged [AP04 eq.n]; everything else is derived here and checked numerically."""
import numpy as np
from scipy.optimize import brentq

# ---- baseline (thesis Table 4.2) ----
P0 = dict(h=145.0, T0=22.0, R=167.0, beta=5.6, eL=0.04, m=0.0, EA=51700.0, EM=46000.0, nu=0.3)

# ---- tensor-level check of the limit function F = sqrt(2 J2) + m J3/J2 - R   [AP04 eq.14] ----
def F_of_X(X, m, R):
    J2 = 0.5*np.trace(X@X); J3 = np.trace(X@X@X)/3.0
    return np.sqrt(2*J2) + m*J3/J2 - R

ntens = np.diag([1.0, -0.5, -0.5]); ntens /= np.linalg.norm(ntens)      # unit tension-like direction (traceless)
kappa_num = (np.trace(ntens@ntens@ntens)/3.0)/(0.5)                       # J3/J2 of a unit tension-like tensor
print("kappa (J3/J2 of unit tension-like tensor) =", kappa_num, " analytic sqrt(2/27)*... =", np.sqrt(2/27)*np.sqrt(3)/np.sqrt(3))

def thresh(sign, m, R):
    """norm b of X = sign*b*ntens at which F=0"""
    return brentq(lambda b: F_of_X(sign*b*ntens, m, R), 1e-6, 10*R)

for m in (0.0, 0.05, 0.459):
    bf = thresh(+1, m, P0['R']); br = thresh(-1, m, P0['R'])
    print(f"m={m}: forward (tension-like X) b = {bf:.3f} = R/(1+k m) -> {P0['R']/(1+kappa_num*m):.3f};  reverse (X=-b n) b = {br:.3f} = R/(1-k m) -> {P0['R']/(1-kappa_num*m):.3f}")
    # uniaxial tension / compression thresholds from the actual stress tensor
    def uni(sig):
        s = sig*np.diag([2/3, -1/3, -1/3]); return F_of_X(s, m, P0['R'])
    st = brentq(uni, 1e-3, 1e4); sc = -brentq(lambda x: uni(-x), 1e-3, 1e4)
    print(f"      uniaxial onset tension {st:.2f} MPa, compression {abs(sc):.2f} MPa;  closed forms {P0['R']/(np.sqrt(2/3)+2*m/9):.2f}, {P0['R']/(np.sqrt(2/3)-2*m/9):.2f}")
    if m > 0:
        sc_, st_ = abs(sc), st
        print(f"      AP04 eq.15 check: R = {2*np.sqrt(2/3)*sc_*st_/(sc_+st_):.2f}, m = {np.sqrt(27/2)*(sc_-st_)/(sc_+st_):.4f}")

# paper's own example set (p.819): sigma_c=72, sigma_t=56
sc_, st_ = 72., 56.
Rp = 2*np.sqrt(2/3)*sc_*st_/(sc_+st_); mp = np.sqrt(27/2)*(sc_-st_)/(sc_+st_)
print(f"\nAP04 p.819 example: R = {Rp:.2f} MPa, m = {mp:.3f} (limit 0.46)")
hp, bp, Mfp, eLp = 1000., 2.1, 223., 0.04
fac = 1 - kappa_num*mp
Af_p = Mfp + Rp/(bp*fac); As_p = Mfp + (Rp/fac - hp*eLp)/bp
print(f"   model-implied zero-stress As = {As_p:.1f} K, Af = {Af_p:.1f} K, window {Af_p-As_p:.1f} K  (listed in paper: As=248, Af=260 K -> inputs for experiment, not outputs of these formulas)")

# ---- free-wire transformation temperatures, window ----
def free(p):
    k = kappa_num
    Af = p['T0'] + p['R']/(p['beta']*(1-k*p['m']))
    As = p['T0'] + (p['R']/(1-k*p['m']) - p['h']*p['eL'])/p['beta']
    return As, Af, Af-As
As0, Af0, W0 = free(P0)
print(f"\nBaseline: As={As0:.3f}, Af={Af0:.3f}, W={W0:.3f}")
rows = []
for name, unit in [('h','MPa'),('beta','MPa/K'),('R','MPa'),('T0','C'),('eL','-'),('m','-')]:
    base = P0[name] if P0[name] != 0 else 0.05
    out = []
    for f in (-0.2, +0.2):
        p = dict(P0)
        if name == 'm': p['m'] = 0.05*(1+f)*1.0 if False else (0.05 if f > 0 else 0.0)
        else: p[name] = P0[name]*(1+f)
        out.append(free(p))
    rows.append((name, out))
    print(f"{name:5s} -20%: dAs={out[0][0]-As0:+7.2f} dAf={out[0][1]-Af0:+7.2f} dW={out[0][2]-W0:+7.3f} | +20%: dAs={out[1][0]-As0:+7.2f} dAf={out[1][1]-Af0:+7.2f} dW={out[1][2]-W0:+7.3f}")
# elastic constants: zero effect (loads s=0)
print("EA, EM, nu: As/Af/W unchanged (stress-free: s = 0 so X = -back stress, moduli never enter)")
# local derivatives
b, h, R, eL, T0, k = P0['beta'], P0['h'], P0['R'], P0['eL'], P0['T0'], kappa_num
print("\nLocal derivatives at baseline (m=0):")
print(f" d(As)/dR=d(Af)/dR = {1/b:.4f} K/MPa;  d(As)/dT0 = d(Af)/dT0 = 1")
print(f" dAs/dbeta = {-(R-h*eL)/b**2:.3f}, dAf/dbeta = {-R/b**2:.3f}, dW/dbeta = {-h*eL/b**2:.3f}  K per (MPa/K)")
print(f" dAs/dh = {-eL/b:.5f}, dW/dh = {eL/b:.5f} K/MPa;  dAs/deL = {-h/b:.2f}, dW/deL = {h/b:.2f} K per unit")
print(f" dAs/dm = dAf/dm = {R*k/b:.2f} K per unit m")

# ---- loaded wire: shift of transformation temperatures with uniaxial stress (m=0) ----
print(f"\nStress shift: dT/dsigma = sqrt(2/3)/beta = {np.sqrt(2/3)/b:.4f} K/MPa = 1 K per {np.sqrt(1.5)*b:.2f} MPa")

# ---- uniaxial stress-strain features (m=0), T = T0 ----
print("\nUniaxial loading at T=T0 (m=0):")
print(f" elastic slope = E_A = {P0['EA']:.0f} MPa; onset = {np.sqrt(1.5)*R:.1f} MPa; plateau rise = {np.sqrt(1.5)*h*eL:.2f} MPa; finish = {np.sqrt(1.5)*(R+h*eL):.1f} MPa;"
      f" transformation strain (uniaxial) = {np.sqrt(2/3)*eL*100:.2f} %; unloading slope after full transformation = E_M (if D(xi) mixing) = {P0['EM']:.0f} MPa")
print(f" plateau tangent slope in total strain = {np.sqrt(1.5)*h*eL/(np.sqrt(2/3)*eL):.1f} MPa (= 1.5 h)  [elastic part not included]")

# ---- constrained (embedded) wire: 1-D reduction, m=0  [derived here] ----
def constrained(p, k_eq, n=2001):
    xi = np.linspace(1, 0, n)
    E = p['EA'] + xi*(p['EM']-p['EA'])
    eps_tr = np.sqrt(2/3)*p['eL']*xi
    sig = (np.sqrt(2/3)*p['eL'] - eps_tr)/(1/E + 1/k_eq)
    T = p['T0'] + (p['R'] + np.sqrt(2/3)*sig - p['h']*p['eL']*xi)/p['beta']
    return xi, sig, T
# calibrate k_eq from supervisor's Af_sigma = 69 C  (As_sigma ~ 51 C is the free value)
Af_sup = 69.0
sig_end = (Af_sup - P0['T0'])*P0['beta'] - P0['R']; sig_end /= np.sqrt(2/3)
c = np.sqrt(2/3)*P0['eL']
Eeff = sig_end/c; k_cal = P0['EA']*Eeff/(P0['EA']-Eeff)
print(f"\nSupervisor Af_sigma=69 C -> end stress {sig_end:.1f} MPa -> series stiffness k_eq = {k_cal:.0f} MPa (referred to wire strain)")
xi, sig, T = constrained(P0, k_cal)
print(f" check: As_sigma = {T[0]:.2f}, Af_sigma = {T[-1]:.2f}, window = {T[-1]-T[0]:.2f} K, end stress = {sig[-1]:.1f} MPa")
print("\nConstrained window vs elastic constants (k_eq fixed):")
for lab, kw in [("baseline", {}), ("EA +20%", dict(EA=P0['EA']*1.2)), ("EA -20%", dict(EA=P0['EA']*0.8)),
                ("EM +20%", dict(EM=P0['EM']*1.2)), ("EM -20%", dict(EM=P0['EM']*0.8)), ("EM = EA", dict(EM=P0['EA']))]:
    p = dict(P0); p.update(kw); xi_, s_, T_ = constrained(p, k_cal)
    mid = np.argmin(abs(xi_-0.5))
    print(f" {lab:9s}: Af_sigma={T_[-1]:.2f} C  window={T_[-1]-T_[0]:.2f} K  end stress={s_[-1]:6.1f} MPa  T(xi=0.5)={T_[mid]:.2f} C  sigma(xi=0.5)={s_[mid]:.1f} MPa")
print("Constrained window vs k_eq (EA,EM baseline):")
for kk in (500, 1000, 2000, 3877, 8000, 20000, 1e9):
    xi_, s_, T_ = constrained(P0, kk)
    print(f" k_eq={kk:>9.0f} MPa: window={T_[-1]-T_[0]:7.2f} K, end stress={s_[-1]:8.1f} MPa")
print(f" dWindow/dEA at k_cal: {( (lambda a,b_: (a[2][-1]-a[2][0])-(b_[2][-1]-b_[2][0]))(constrained(dict(P0,EA=P0['EA']*1.01),k_cal), constrained(P0,k_cal)) )/0.01:.2f} K per unit relative change")

# ---- Clausius-Clapeyron / entropy consistency of beta ---- [AP04 eq.8_3 : entropy term -beta*||e_tr||]
ds = b*eL                      # MPa/K * 1 = J/mm^3/K *1e-3 -> mJ/mm3/K
for Tt in (As0+273.15, Af0+273.15):
    L = Tt*ds*1e-3             # J/mm3
    print(f"\nModel-implied latent heat at T={Tt:.1f} K: {L:.4f} J/mm3 = {L*1e9/6450/1e3:.1f} kJ/kg ; wire V=628 mm3 -> {L*628:.1f} J (plant uses 0.156 J/mm3 = {0.156*628:.1f} J)")
