"""Brute-force search for correct z0^CM at p=11, s=1/2
Test various candidate forms:
- 1/p^k for k=2..8
- alpha/p^2 for rational alpha in [0.5, 2]
- Known Heegner j-invariant based formulas
"""
from sage.all import *
RR = RealField(500)
pi_val = RR.pi()

a = [RR(1)/2, RR(1)/22, RR(3)/22, RR(5)/22, RR(7)/22, RR(9)/22,
     RR(13)/22, RR(15)/22, RR(17)/22, RR(19)/22, RR(21)/22]
b = [RR(1)/11, RR(2)/11, RR(3)/11, RR(4)/11, RR(5)/11,
     RR(6)/11, RR(7)/11, RR(8)/11, RR(9)/11, RR(10)/11]

def poch(r, n):
    return gamma(r + n) / gamma(r)

def compute_S(z0, max_n=50):
    S0 = RR(0); S1 = RR(0)
    for n in range(max_n):
        num = prod(poch(ai, n) for ai in a)
        den = prod(poch(bj, n) for bj in b)
        den *= gamma(RR(n + 1))
        r = num / den
        S0 += r * z0**n
        S1 += r * RR(n) * z0**n
        if n > 2 and abs(r*z0**n)/abs(S0) < RR(2)**(-500):
            break
    return S0, S1

def scan_z0(z0, label, max_A=100, max_B=1000):
    S0, S1 = compute_S(z0)
    pS0 = pi_val * S0
    pS1 = pi_val * S1
    hits = []
    for A in range(1, max_A+1):
        for B in range(1, max_B+1):
            C_sq = (RR(A)*pS0 + RR(B)*pS1)**2
            rat = C_sq.nearby_rational(max_denominator=100000)
            err = abs(C_sq - RR(rat))
            if err < RR('1e-20'):
                hits.append((float(err), A, B, rat, float(pS0), float(pS1)))
    return hits, pS0, pS1

# Candidates based on various formulas
candidates = {}
# Power of 11
for k in range(2, 9):
    candidates[f"1/11^{k}"] = RR(1) / (RR(11)**k)
# Rational multiples of 1/121
for num in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]:
    candidates[f"{num}/121"] = RR(num) / 121
# j(i)=1728 based formula from paper: 1728^{-1/11}/121
candidates["1728^{-1/11}/121"] = RR(1728)**(RR(-1)/11) / RR(121)
# lambda(11i)
candidates["lambda(11i)"] = RR(1.5702908149651593e-14)
# 1/11^2 * sqrt(3) etc
candidates["sqrt(3)/121"] = RR(3).sqrt() / 121
# (sqrt(11)-1)/121
candidates["(sqrt(11)-3)/121"] = (RR(11).sqrt() - 3) / 121
candidates["(sqrt(11)-1)/121"] = (RR(11).sqrt() - 1) / 121

print(f"{'Candidate':<30} {'A':<4} {'B':<5} {'C^2':<25} {'Error':<12}")
print(f"{'-'*30} {'-'*4} {'-'*5} {'-'*25} {'-'*12}")
for lbl, z0 in candidates.items():
    hits, _, _ = scan_z0(z0, lbl, max_A=50, max_B=500)
    if hits:
        hits.sort(key=lambda x: x[0])
        for err, A, B, C2rat, _, _ in hits[:3]:
            print(f"{lbl:<30} {A:<4} {B:<5} {str(C2rat):<25} {err:.2e}")
    else:
        print(f"{lbl:<30} {'--':<4} {'--':<5} {'NO HITS':<25}")

# Also try a fine-grained search near 1/121
print(f"\nFine scan: z0 = alpha/121 for alpha in [0.1, 2.0] step 0.05")
for alpha_int in range(10, 201, 5):
    alpha = RR(alpha_int) / 100
    z0 = alpha / 121
    hits, pS0, pS1 = scan_z0(z0, f"alpha={alpha}", max_A=30, max_B=200)
    if hits:
        hits.sort(key=lambda x: x[0])
        err, A, B, C2rat, _, _ = hits[0]
        if err < RR('1e-25'):
            print(f"  alpha={alpha:.2f}: (A,B)=({A},{B}), C^2={C2rat}, error={err:.2e}")

print("\nDone.")
