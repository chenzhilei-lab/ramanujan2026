"""Verify suspicious PSLQ hits from first scan with tighter threshold"""
from sage.all import *
RR = RealField(500)
pi_val = RR.pi()

a = [RR(1)/2, RR(1)/22, RR(3)/22, RR(5)/22, RR(7)/22, RR(9)/22,
     RR(13)/22, RR(15)/22, RR(17)/22, RR(19)/22, RR(21)/22]
b = [RR(1)/11, RR(2)/11, RR(3)/11, RR(4)/11, RR(5)/11,
     RR(6)/11, RR(7)/11, RR(8)/11, RR(9)/11, RR(10)/11]

def poch(r, n):
    return gamma(r + n) / gamma(r)

# Suspicious z0 values from first scan
candidates = [
    (43, 0.99928534992, 20, 441),
    (43, 0.5, 48, 471),
    (43, 0.36608682733, 8, 282),
    (67, 0.5, 48, 471),
    (67, 0.89992577069, 5, 120),
    (67, 0.22783842379, 8, 458),
]

print(f"{'d':<5} {'z0':<20} {'A':<5} {'B':<5} {'C^2_err (500-bit)':<20}")
print(f"{'-'*5} {'-'*20} {'-'*5} {'-'*5} {'-'*20}")

for d, z0_val, A, B in candidates:
    z0 = RR(z0_val)
    S0 = RR(0); S1 = RR(0)
    for n in range(100):
        num = prod(poch(ai, n) for ai in a)
        den = prod(poch(bj, n) for bj in b)
        den *= gamma(RR(n + 1))
        r = num / den
        S0 += r * z0**n
        S1 += r * RR(n) * z0**n
        if n > 2 and abs(r*z0**n)/abs(S0) < RR(2)**(-400):
            break
    
    C_sq = (RR(A)*pi_val*S0 + RR(B)*pi_val*S1)**2
    rat = C_sq.nearby_rational(max_denominator=100000)
    err = abs(C_sq - RR(rat))
    if err < RR('1e-30'):
        print(f"{d:<5} {z0_val:<20.10e} {A:<5} {B:<5} {err:<20.2e}  GENUINE?")
    else:
        print(f"{d:<5} {z0_val:<20.10e} {A:<5} {B:<5} {err:<20.2e}  NOISE (>{'1e-30':>6})")

# Also compute exact pi*S0, pi*S1 for each z0 and check if the C^2 ratio is 
# a simple algebraic number
print(f"\n{'='*70}")
print("Detailed check: does C^2 match a rational with small denominator?")
print(f"{'='*70}")
print(f"{'z0':<20} {'(A,B)':<12} {'C^2_as_rational':<30} {'denom':<10}")
for d, z0_val, A, B in candidates:
    z0 = RR(z0_val)
    S0 = RR(0); S1 = RR(0)
    for n in range(100):
        num = prod(poch(ai, n) for ai in a)
        den = prod(poch(bj, n) for bj in b)
        den *= gamma(RR(n + 1))
        r = num / den
        S0 += r * z0**n
        S1 += r * RR(n) * z0**n
        if n > 2 and abs(r*z0**n)/abs(S0) < RR(2)**(-400):
            break
    C_sq = (RR(A)*pi_val*S0 + RR(B)*pi_val*S1)**2
    rat = C_sq.nearby_rational(max_denominator=100000)
    err = abs(C_sq - RR(rat))
    print(f"{z0_val:<20.10e} ({A},{B})    {str(rat):<30} {rat.denom():<10} err={err:.2e}")

print("\nDone.")
