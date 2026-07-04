"""PSLQ search at candidate singular moduli for F1"""
from sage.all import *

RR = RealField(500)
pi_val = RR.pi()

# F1 parameters
a = [RR(1)/2, RR(1)/22, RR(3)/22, RR(5)/22, RR(7)/22, RR(9)/22,
     RR(13)/22, RR(15)/22, RR(17)/22, RR(19)/22, RR(21)/22]
b = [RR(1)/11, RR(2)/11, RR(3)/11, RR(4)/11, RR(5)/11,
     RR(6)/11, RR(7)/11, RR(8)/11, RR(9)/11, RR(10)/11]

def poch(r, n):
    return gamma(r + n) / gamma(r)

def compute_S(z0, max_n=100):
    S0 = RR(0); S1 = RR(0)
    for n in range(max_n):
        num = prod(poch(ai, n) for ai in a)
        den = prod(poch(bj, n) for bj in b)
        den *= gamma(RR(n + 1))
        r = num / den
        term0 = r * z0**n
        term1 = r * RR(n) * z0**n
        S0 += term0; S1 += term1
        if n > 2 and abs(term0)/abs(S0) < RR(2)**(-500):
            break
    return S0, S1

def pslq_scan(z0, name, max_A=50, max_B=500):
    S0, S1 = compute_S(z0)
    piS0 = pi_val * S0
    piS1 = pi_val * S1
    
    hits = []
    for A in range(1, max_A+1):
        for B in range(1, max_B+1):
            C_sq = (RR(A)*piS0 + RR(B)*piS1)**2
            rat = C_sq.nearby_rational(max_denominator=100000)
            err = abs(C_sq - RR(rat))
            if err < RR('1e-30'):
                hits.append((float(err), A, B, rat, float(piS0), float(piS1)))
    
    return hits, piS0, piS1

# Candidate moduli to test
trials = {
    "1/11^2 = 1/121": RR(1)/121,
    "1/11^3 = 1/1331": RR(1)/1331,
    "1/11^4 = 1/14641": RR(1)/14641,
    "lambda(11i)": RR(1.5702908149651593e-14),
}

print(f"{'Candidate':<25} {'A':<4} {'B':<5} {'C^2':<25} {'Error':<12}")
print(f"{'-'*25} {'-'*4} {'-'*5} {'-'*25} {'-'*12}")

for name, z0_val in trials.items():
    hits, piS0, piS1 = pslq_scan(z0_val, name)
    if hits:
        hits.sort(key=lambda x: x[0])
        err, A, B, C2rat, _, _ = hits[0]
        print(f"{name:<25} {A:<4} {B:<5} {str(C2rat):<25} {err:.2e}")
    else:
        print(f"{name:<25} {'--':<4} {'--':<5} {'NO HITS':<25} {'':<12}")

# Also try a numerical search with wider range at 1/121
print(f"\n\nDetailed scan at z0 = 1/121 (trial modulus):")
hits121, piS0_121, piS1_121 = pslq_scan(RR(1)/121, max_A=100, max_B=1000)
if hits121:
    hits121.sort(key=lambda x: x[0])
    print(f"{'A':<4} {'B':<5} {'C^2':<30} {'Error':<12}")
    for i, (err, A, B, C2rat, _, _) in enumerate(hits121[:10]):
        print(f"{A:<4} {B:<5} {str(C2rat):<30} {err:.2e}")
else:
    print("No hits at 1/121 either!")

print(f"\npi*S0 = {piS0_121}")
print(f"pi*S1 = {piS1_121}")
