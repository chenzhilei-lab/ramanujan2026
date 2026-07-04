"""
Quick scan: find z0 that gives clean S0*pi for F1 (s=1/2, p=11).
Pattern from known s=1/2: z0 = 1/p^p (p=2→1/4, p=3→1/27, p=5→1/3125).
"""
from mpmath import mp, mpf, sqrt, pi, gamma
from fractions import Fraction

mp.dps = 200

a = [Fraction(1,2), Fraction(1,22), Fraction(3,22), Fraction(5,22),
     Fraction(7,22), Fraction(9,22), Fraction(13,22), Fraction(15,22),
     Fraction(17,22), Fraction(19,22), Fraction(21,22)]
b = [Fraction(1,11), Fraction(2,11), Fraction(3,11), Fraction(4,11),
     Fraction(5,11), Fraction(6,11), Fraction(7,11), Fraction(8,11),
     Fraction(9,11), Fraction(10,11)]

def poch(r, n):
    x = mpf(r.numerator) / r.denominator
    return gamma(x + n) / gamma(x)

def R(n, z):
    num = mpf(1)
    for ai in a: num *= poch(ai, n)
    den = mpf(1)
    for bj in b: den *= poch(bj, n)
    for _ in range(len(a) - len(b)): den *= gamma(n + 1)
    return num / den * z**n

def analyze(z0, label):
    S0, S1 = mpf(0), mpf(0)
    for n in range(200):
        rn = R(n, z0)
        S0_old, S1_old = S0, S1
        S0 += rn; S1 += rn * n
        if n > 5 and abs(rn) / max(abs(S0), mpf('1e-80')) < mpf('1e-60'):
            break
    piS0 = S0 * pi
    print(f"  {label:20s}  S0={float(S0):.12f}  S0*pi={float(piS0):.12f}  "
          f"terms={n}  S1/S0={float(S1/S0):.6e}")

print("=" * 65)
print("F1 (s=1/2, p=11): z0 scan")
print("=" * 65)

# Pattern: z0 = 1/p^p
analyze(mpf(1)/mpf(11)**11, 'z0=1/11^11')
analyze(mpf(1)/mpf(11)**2, 'z0=1/11^2')
analyze(mpf(1)/mpf(11)**3, 'z0=1/11^3')
analyze(mpf(1)/mpf(11)**4, 'z0=1/11^4')
analyze(mpf(1)/mpf(11)**5, 'z0=1/11^5')
analyze(mpf(1)/mpf(11)**6, 'z0=1/11^6')

# Also try based on possible CM discriminant scaling
# For d-based: z0 = 1/(small integer)^{exponent}
for d in [1, 2, 3, 7, 11, 19]:
    j_vals = {1: 1728, 2: 8000, 3: 0, 7: -3375, 11: -32768, 19: -884736}
    jv = j_vals[d]
    if jv == 0: continue
    # The correct z0 at d might be related to the singular modulus
    # For d=1: try z0 = 1/4 (known for p=2) and scaled versions

print(f"\n  For d=1 (tau=i): trying Ramanujan-style z0 candidates")

# The singular modulus for level 2 at tau=i: lambda(i) = 1/2
# Possibly z0 is 1/(some function of 1/2 and p)
# For p=2: z0 = 1/4 = lambda(i)^2 = (1/2)^2
# For p=3: z0 = 1/27 — not a simple power of 1/2
# The level-p modular equation changes the value

# Try the pattern: z0 = (lambda(i))^p / p^2 = (1/2)^p / p^2
for p_try in [11, 7, 5, 3]:
    z_try = mpf(0.5)**p_try / mpf(p_try*p_try)
    analyze(z_try, f'(1/2)^{p_try}/{p_try}^2')

# Try: z0 = 1 / (CM_discriminant_absolute_j)^{exponent}
# For |j| = 1728 at d=1:
# If z0 = 1 / 1728^K, find K
for K in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
    z_try = mpf(1) / (mpf(1728)**K)
    analyze(z_try, f'1/1728^{K}')

print(f"\n  Current best: 1/11^6 gives S0={float(mpf(1)/mpf(11)**6):.6e} z0")
