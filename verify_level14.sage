"""Verify the level-14 formula claimed in the framework paper.
s=1/2, p=7, d=7, z0=1/196, k=7, 7 numerator params
"""
from sage.all import *
RR = RealField(500)
pi_val = RR.pi()

# Level-14 numerator parameters (from framework paper)
a = [RR(1)/14, RR(3)/14, RR(5)/14, RR(7)/14, RR(9)/14, RR(11)/14, RR(13)/14]
# Denominator: (n!)^7 means |a'|-|b'| = k = 7, so |b'| = 0
b = []  # No explicit denominator parameters; all in (n!)^7

z0 = RR(1)/196

def poch(r, n):
    return gamma(r + n) / gamma(r)

print("Level-14 formula verification (s=1/2, p=7, d=7):")
print(f"  |a'| = {len(a)}, k = 7, z0 = 1/196")
print(f"  Series: sum prod(a_i)_n / (n!)^7 * (A+Bn) * (1/196)^n = C/pi")
print()

S0 = RR(0); S1 = RR(0)
for n in range(100):
    num = prod(poch(ai, n) for ai in a)
    den = gamma(RR(n+1))**7
    r = num / den
    S0 += r * z0**n
    S1 += r * RR(n) * z0**n
    if n > 5 and abs(r*z0**n)/abs(S0) < RR(2)**(-400):
        print(f"  Converged at n={n+1}")
        break

pS0 = pi_val * S0
pS1 = pi_val * S1
print(f"\n  pi*S0 = {pS0}")
print(f"  pi*S1 = {pS1}")

# PSLQ search A in [1,100], B in [1,500]
print(f"\n  PSLQ search A=[1,100], B=[1,500]...")
hits = []
for A in range(1, 101):
    for B in range(1, 501):
        C_sq = (RR(A)*pS0 + RR(B)*pS1)**2
        rat = C_sq.nearby_rational(max_denominator=100000)
        err = abs(C_sq - RR(rat))
        if err < RR('1e-30'):
            hits.append((float(err), A, B, rat))

if hits:
    hits.sort(key=lambda x: x[0])
    print(f"  Found {len(hits)} hits!")
    print(f"  {'A':<5} {'B':<5} {'C^2':<30} {'Error':<12}")
    for err, A, B, C2rat in hits[:10]:
        print(f"  {A:<5} {B:<5} {str(C2rat):<30} {err:.2e}")
else:
    print("  NO HITS — formula does NOT exist at these parameters")
    
    # Try wider search
    print(f"\n  Extended search A=[1,200], B=[1,2000]...")
    for A in range(1, 201):
        for B in range(1, 2001):
            C_sq = (RR(A)*pS0 + RR(B)*pS1)**2
            rat = C_sq.nearby_rational(max_denominator=200000)
            err = abs(C_sq - RR(rat))
            if err < RR('1e-30'):
                hits.append((float(err), A, B, rat))

    if hits:
        hits.sort(key=lambda x: x[0])
        print(f"  Found {len(hits)} hits in extended range!")
        for err, A, B, C2rat in hits[:10]:
            print(f"  {A:<5} {B:<5} {str(C2rat):<30} {err:.2e}")
    else:
        print("  STILL NO HITS — level-14 formula is NOT valid")

print("\nDone.")
