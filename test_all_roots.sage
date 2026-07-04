"""Compute ALL roots of the modular equation and PSLQ test each lambda value"""
from sage.all import *

# Get modular polynomial for j-function at level 11
phi = pari.polmodular(11)

# Build Phi_11(x, 1728) = 0 polynomial
poly_x = pari(0)
var_x = pari('x')
var_y = pari('y')
for k in range(13):
    coeff = pari.subst(phi[k], var_y, 1728)
    poly_x = poly_x + coeff * var_x**k

# Find ALL roots 
roots = pari.polroots(poly_x)
print(f"Total roots of Phi_11(x, 1728) = 0: {len(roots)}")

# Convert each j-value to lambda and test with PSLQ
RR = RealField(500)
CC = ComplexField(500)
pi_val = RR.pi()

# F1 parameters
a = [RR(1)/2, RR(1)/22, RR(3)/22, RR(5)/22, RR(7)/22, RR(9)/22,
     RR(13)/22, RR(15)/22, RR(17)/22, RR(19)/22, RR(21)/22]
b = [RR(1)/11, RR(2)/11, RR(3)/11, RR(4)/11, RR(5)/11,
     RR(6)/11, RR(7)/11, RR(8)/11, RR(9)/11, RR(10)/11]

def poch(r, n):
    return gamma(r + n) / gamma(r)

def test_z0(z0_val, max_A=50, max_B=500):
    S0 = RR(0); S1 = RR(0)
    for n in range(100):
        num = prod(poch(ai, n) for ai in a)
        den = prod(poch(bj, n) for bj in b)
        den *= gamma(RR(n + 1))
        r = num / den
        S0 += r * z0_val**n
        S1 += r * RR(n) * z0_val**n
        if n > 2 and abs(r*z0_val**n)/abs(S0) < RR(2)**(-500):
            break
    
    pS0 = pi_val * S0
    pS1 = pi_val * S1
    
    best = None
    for A in range(1, max_A+1):
        for B in range(1, max_B+1):
            C_sq = (RR(A)*pS0 + RR(B)*pS1)**2
            rat = C_sq.nearby_rational(max_denominator=100000)
            err = abs(C_sq - RR(rat))
            if err < RR('1e-20'):
                if best is None or err < best[0]:
                    best = (float(err), A, B, rat)
    return best

# Test each root
print(f"\n{'Root #':<7} {'|j|':<15} {'arg(j)':<12} {'lambda':<15} {'Best (A,B)':<15} {'Error':<12}")
print(f"{'-'*7} {'-'*15} {'-'*12} {'-'*15} {'-'*15} {'-'*12}")

for i in range(12):
    j_val = roots[i]
    j_real = float(pari.real(j_val))
    j_imag = float(pari.imag(j_val))
    j_abs = sqrt(j_real**2 + j_imag**2)
    j_arg = float(pari.arg(j_val)) if hasattr(pari, 'arg') else 0
    
    # Convert j to lambda using Newton's method
    # j = 256*(1-lambda+lambda^2)^3/(lambda^2*(1-lambda)^2)
    # Solve for lambda using inverse formula
    
    # For small |j| use Newton from lam=0.5
    # For large |j| use lam = 16/sqrt(j) as initial (small lambda approximation)
    
    try:
        if j_abs < 1e10:
            lam_guess = RR(0.5)
        else:
            C = CC(j_val)
            lam_guess = 16 / C.sqrt()
        
        # Newton refinement
        lam = CC(lam_guess)
        for _ in range(30):
            l = lam
            c_j = CC(j_val)
            f = c_j * l*l * (1-l)*(1-l) - 256 * (1 - l + l*l)**3
            df = c_j * (2*l*(1-l)**2 - 2*l*l*(1-l)) - 768 * (1-l+l*l)**2 * (-1+2*l)
            if abs(df) < 1e-100: break
            lam = l - f/df
        
        lam_real = abs(lam.real())
        
        # Only test if lambda is roughly in [0,1]
        if 0 < lam_real < 1:
            best = test_z0(RR(lam_real))
        else:
            best = None
    except:
        best = None
    
    if best:
        err, A, B, C2rat = best
        print(f"{j_abs:<15.6e} {j_arg:<12.4f} {lam_real:<15.10e} ({A},{B})  {err:.2e}")
    elif 0 < (lam_real if 'lam_real' in dir() else -1) < 1:
        print(f"{j_abs:<15.6e} {j_arg:<12.4f} {lam_real:<15.10e} {'--':>15}  NO HITS")
    else:
        pass  # Skip unstable roots

print("\nDone.")
