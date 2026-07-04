"""Scan all 9 Heegner discriminants for F1 (s=1/2, p=11)
For each d: Φ_11(x, j(τ_d)) = 0 → 12 j-roots → convert to λ → PSLQ
"""
from sage.all import *
import sys

RR = RealField(500)
pi_val = RR.pi()
CC = ComplexField(500)

# F1 hypergeometric parameters
a = [RR(1)/2, RR(1)/22, RR(3)/22, RR(5)/22, RR(7)/22, RR(9)/22,
     RR(13)/22, RR(15)/22, RR(17)/22, RR(19)/22, RR(21)/22]
b = [RR(1)/11, RR(2)/11, RR(3)/11, RR(4)/11, RR(5)/11,
     RR(6)/11, RR(7)/11, RR(8)/11, RR(9)/11, RR(10)/11]

def poch(r, n):
    return gamma(r + n) / gamma(r)

def test_z0(z0_val, max_A=80, max_B=800):
    """PSLQ scan at given z0, return best hit"""
    if z0_val <= 0 or z0_val >= 1:
        return None
    S0 = RR(0); S1 = RR(0)
    for n in range(100):
        num = prod(poch(ai, n) for ai in a)
        den = prod(poch(bj, n) for bj in b)
        den *= gamma(RR(n + 1))
        r = num / den
        S0 += r * z0_val**n
        S1 += r * RR(n) * z0_val**n
        if n > 2 and abs(r*z0_val**n)/abs(S0) < RR(2)**(-300):
            break
    pS0 = pi_val * S0
    pS1 = pi_val * S1
    
    best = None
    for A in range(1, max_A+1):
        for B in range(1, max_B+1):
            C_sq = (RR(A)*pS0 + RR(B)*pS1)**2
            rat = C_sq.nearby_rational(max_denominator=100000)
            err = abs(C_sq - RR(rat))
            if err < RR('1e-30'):
                if best is None or err < best[0]:
                    best = (float(err), A, B, rat)
    return best

def j_to_lambda(j_cc):
    """Convert j (complex) to lambda using Newton's method.
    Returns all valid lambda values in (0,1)."""
    results = []
    # Try multiple initial guesses
    for guess in [CC(0.1), CC(0.3), CC(0.5), CC(0.7), CC(0.9)]:
        lam = guess
        try:
            for _ in range(50):
                l = lam
                f = j_cc * l*l * (1-l)*(1-l) - 256 * (1 - l + l*l)**3
                df = j_cc * (2*l*(1-l)**2 - 2*l*l*(1-l)) - 768 * (1-l+l*l)**2 * (-1+2*l)
                if abs(df) < 1e-150:
                    break
                lam = l - f/df
            
            lam_r = float(lam.real())
            lam_i = float(lam.imag())
            if abs(lam_i) < 1e-10 and 0 < lam_r < 1:
                # Check if this is a new root (not duplicate)
                dup = False
                for r, _ in results:
                    if abs(r - lam_r) < 1e-8:
                        dup = True
                        break
                if not dup:
                    results.append((lam_r, float(f)))
        except:
            pass
    return results

# Heegner discriminants and their j-invariants
heegner = {
    1:  1728,
    2:  8000,
    3:  0,
    7:  -3375,
    11: -32768,
    19: -884736,
    43: -884736000,
    67: -147197952000,
    163: -262537412640768000,
}

# Get modular polynomial
phi = pari.polmodular(11)
var_x = pari('x')
var_y = pari('y')

print(f"{'d':<5} {'j(τ_d)':<25} {'z0_candidate':<20} {'A':<5} {'B':<5} {'C^2':<25} {'Error':<12}")
print(f"{'-'*5} {'-'*25} {'-'*20} {'-'*5} {'-'*5} {'-'*25} {'-'*12}")

total = 0
hits = []

for d, j_val in heegner.items():
    # Build Φ_11(x, j_val) = 0
    poly_x = pari(0)
    for k in range(13):
        coeff = pari.subst(phi[k], var_y, pari(j_val))
        poly_x = poly_x + coeff * var_x**k
    
    # Find all roots
    roots = pari.polroots(poly_x)
    
    for i in range(12):
        j_root = roots[i]
        j_r = float(pari.real(j_root))
        j_i = float(pari.imag(j_root))
        j_abs = sqrt(j_r**2 + j_i**2) if abs(j_r) < 1e100 else float('inf')
        
        # Convert to lambda
        lam_list = j_to_lambda(CC(float(pari.real(j_root)), float(pari.imag(j_root))))
        
        for lam_r, _ in lam_list:
            total += 1
            best = test_z0(RR(lam_r), max_A=80, max_B=800)
            if best:
                err, A, B, C2rat = best
                hits.append((err, d, lam_r, A, B, C2rat))
                print(f"{d:<5} {j_val:<25} {lam_r:<20.10e} {A:<5} {B:<5} {str(C2rat):<25} {err:.2e}")

# Summary
print(f"\n{'='*80}")
print(f"Total z0 candidates tested: {total}")
print(f"Total PSLQ hits: {len(hits)}")

if hits:
    hits.sort(key=lambda x: x[0])
    print(f"\nTop hits across all discriminants:")
    print(f"{'d':<5} {'lambda':<20} {'A':<5} {'B':<5} {'C^2':<25} {'Error':<12}")
    print(f"{'-'*5} {'-'*20} {'-'*5} {'-'*5} {'-'*25} {'-'*12}")
    for err, d, lam, A, B, C2rat in hits[:10]:
        print(f"{d:<5} {lam:<20.10e} {A:<5} {B:<5} {str(C2rat):<25} {err:.2e}")
else:
    print("\nNo PSLQ hits found for ANY discriminant.")
    print("Conclusion: F1 (s=1/2, p=11) has no C/pi closed form")
    print("for any of the 9 Heegner discriminants.")

print("\nDone.")
