"""Compute exact CM singular modulus for F1 using PARI's modular polynomial"""
from sage.all import *

# Step 1: Get modular polynomial for j-function at level 11
# polmodular(N) returns a vector of polynomials
# For classical modular polynomial Phi_N(x, y) = 0 where x=j(N*tau), y=j(tau)
phi_data = pari.polmodular(11)
print(f"phi_data length: {len(phi_data)}")

# phi_data is a 13-element vector. The first element is the polynomial in y (j(tau))
# Actually in PARI/GP: polmodular(N) returns a vector [F, G1, G2, ...]
# where F relates j(N*tau) to j(tau)
# Let me check what each element is
print(f"phi_data[0] = {phi_data[0]}")
print(f"poldeg = {pari.poldegree(phi_data[0])}")

# This polynomial Phi_11(y, j) = 0 where y = j(11*tau), j = j(tau)
# Substitute j = j(i) = 1728
R = PolynomialRing(QQ, 'y')
y = R.gen()
Phi11_y = R(pari.subst(phi_data[0], 'x', 1728))
print(f"\nPhi_11(y, 1728) degree: {Phi11_y.degree()}")

# Find real roots
roots_real = []
for root, mult in Phi11_y.roots(RR):
    print(f"  Root: y = {float(root):.6e}")
    roots_real.append(float(root))

# The root j(11i) should be about 1.038e30 based on PARI's computation
target = 1.03819703567651e30
print(f"\nExpected j(11i) ≈ {target:.6e}")

# Find the closest root
closest = min(roots_real, key=lambda r: abs(r - target))
print(f"Closest root: {closest:.6e}")
print(f"Ratio: {closest/target:.6f}")

# Step 2: Convert j to lambda
# j = 256*(1 - lambda + lambda^2)^3 / (lambda^2 * (1-lambda)^2)
# Given j, find lambda by solving j*lambda^2*(1-lambda)^2 = 256*(1-lambda+lambda^2)^3

j11 = closest
R2 = PolynomialRing(QQbar, 'lam')
lam = R2.gen()
eq = j11 * lam**2 * (1-lam)**2 - 256 * (1 - lam + lam**2)**3

# Find roots of this equation
lam_roots = eq.roots(multiplicities=False)
print(f"\nConverting j(11i) to lambda(11i):")
for lr in lam_roots:
    abs_lr = abs(RR(lr))
    print(f"  lambda = {lr}, |lambda| = {float(abs_lr):.6e}")

# The correct lambda should be small (for convergence)
lambdas_small = [(float(abs(lr)), lr) for lr in lam_roots if abs(lr) < 1]
lambdas_small.sort(key=lambda x: x[0])
print(f"\nSmallest lambda (candidate z0^CM): {lambdas_small[0][1]}")
print(f"  |z0^CM| = {lambdas_small[0][0]:.6e}")

z0_exact = lambdas_small[0][1]
print(f"\n  lambda(11i) from theta functions: 1.57029081496516e-14")
print(f"  lambda(11i) from j-polynomial:     {float(abs(z0_exact)):.6e}")

# Step 3: Now compute PSLQ at the exact modulus
print("\n" + "="*70)
print("STEP 3: PSLQ at exact z0^CM")
print("="*70)

R200 = RealField(200)

# F1 parameters
a = [R200(1)/2, R200(1)/22, R200(3)/22, R200(5)/22, R200(7)/22, R200(9)/22,
     R200(13)/22, R200(15)/22, R200(17)/22, R200(19)/22, R200(21)/22]
b = [R200(1)/11, R200(2)/11, R200(3)/11, R200(4)/11, R200(5)/11,
     R200(6)/11, R200(7)/11, R200(8)/11, R200(9)/11, R200(10)/11]

z0 = R200(abs(z0_exact))

def poch(r, n):
    return gamma(r + n) / gamma(r)

S0 = R200(0); S1 = R200(0)
for n in range(100):
    num = prod(poch(ai, n) for ai in a)
    den = prod(poch(bj, n) for bj in b)
    den *= gamma(R200(n + 1))
    r = num / den
    zn = z0**n
    term0 = r * zn
    term1 = r * n * zn
    S0 += term0; S1 += term1
    if n > 3 and abs(term0)/abs(S0) < R200(2)**(-200):
        break

pi_val = R200(pi)
print(f"S0 = {S0}")
print(f"S1 = {S1}")
print(f"pi*S0 = {pi_val*S0}")
print(f"pi*S1 = {pi_val*S1}")

# PSLQ search
print(f"\nSearching A in [1,50], B in [1,500] at exact z0^CM...")
best_err = R200(1); best_AB = (0,0)
hits = []
for A in range(1, 51):
    for B in range(1, 501):
        C_pi = R200(A)*S0 + R200(B)*S1
        C = pi_val * C_pi
        C2 = C*C
        rat = C2.nearby_rational(max_denominator=50000)
        err = abs(C2 - R200(rat))
        if err < R200('1e-30'):
            C_rat = sqrt(rat)
            hits.append((float(err), A, B, rat))
            if float(err) < float(best_err):
                best_err = err; best_AB = (A,B)

if hits:
    hits.sort(key=lambda x: x[0])
    print(f"\nTop hits at exact z0^CM:")
    for i, (err, A, B, C2rat) in enumerate(hits[:5]):
        C_val = sqrt(R200(C2rat))
        print(f"  (A,B)=({A},{B}), C^2={C2rat}, error={err:.2e}")
    
    # Check if (16,250) survives
    C_trial = sqrt(R200(2381505))/R200(30)
    C_pi_check = R200(16)*S0 + R200(250)*S1
    C_check = pi_val * C_pi_check
    print(f"\nTrial candidate (16,250) at exact modulus:")
    print(f"  C = {C_check}")
    print(f"  Expected (if invariant): {C_trial}")
    print(f"  Error: {abs(C_check-C_trial)/abs(C_trial):.6e}")
else:
    print("No PSLQ hits found within search range.")
