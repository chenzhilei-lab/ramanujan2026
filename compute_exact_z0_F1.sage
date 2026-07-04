#!/usr/bin/env sage
"""
Compute exact CM singular modulus for Family F1 (s=1/2, p=11, d=1)
and re-run PSLQ at the exact modulus to find (A*,B*).

Usage:
    sage compute_exact_z0_F1.sage 2>&1 | tee f1_exact_results.txt
"""

# ======================================================================
# Step 1: Get the modular polynomial Phi_11
# ======================================================================
print("=" * 70)
print("STEP 1: Computing modular polynomial Phi_11")
print("=" * 70)

from sage.databases.db_modular_polynomials import ClassicalModularPolynomialDatabase
Phi11 = ClassicalModularPolynomialDatabase()[11]

print(f"Phi_11(X,Y) is a bivariate polynomial of total degree {Phi11.total_degree()}")
print(f"Number of terms: {Phi11.number_of_terms()}")

# ======================================================================
# Step 2: Substitute Y = lambda(i) = 1/2
# ======================================================================
print("\n" + "=" * 70)
print("STEP 2: Substituting Y = lambda(i) = 1/2")
print("=" * 70)

R = QQ['X']
X = R.gen()
phi11_X = R(Phi11(X, QQ(1)/2))

print(f"Degree of Phi_11(X, 1/2): {phi11_X.degree()}")
print(f"Coefficients (top 5): {phi11_X.coefficients()[:5]}...")

# ======================================================================
# Step 3: Find the correct algebraic root
# ======================================================================
print("\n" + "=" * 70)
print("STEP 3: Finding the CM singular modulus")
print("=" * 70)

# Factor the polynomial to find its roots over Qbar
# We need the root that corresponds to z0^CM for the level-raised seed
# For s=1/2, p=11, the correct root is typically the one with 
# |z0| < 1 (for convergence) and with positive real part

# Factor over QQ first
factors = factor(phi11_X)
print(f"Factorization over QQ: {factors}")

# Compute all real roots numerically to find the right one
# (The one with correct convergence properties)
RR = RealField(200)
roots = []
for root, mult in phi11_X.roots(RR):
    roots.append((root, mult))
    print(f"  Real root: root = {root}")

# If there are no real roots, look for complex roots with |z| small
if not roots:
    CC = ComplexField(200)
    for root, mult in phi11_X.roots(CC):
        roots.append((root, mult))
        print(f"  Complex root: root = {root}, |root| = {abs(root):.6e}")

# Sort by absolute value - the singular modulus should be small
roots_sorted = sorted(roots, key=lambda r: abs(r[0]))
print(f"\nSmallest root (candidate z0^CM): {roots_sorted[0][0]}")
print(f"  |z0^CM| = {abs(roots_sorted[0][0]):.6e}")
print(f"  Trial modulus 1/121 = {1/121:.6e}")
print(f"  Ratio |z0^CM| / (1/121) = {abs(roots_sorted[0][0]) / (1/121):.6e}")

z0_exact = roots_sorted[0][0]

# ======================================================================
# Step 4: Compute the series at the exact modulus
# ======================================================================
print("\n" + "=" * 70)
print("STEP 4: Computing series sum at exact z0^CM")
print("=" * 70)

# Use high precision
R200 = RealField(200)

# F1 parameters
a = [R200(1)/2, R200(1)/22, R200(3)/22, R200(5)/22, R200(7)/22, R200(9)/22,
     R200(13)/22, R200(15)/22, R200(17)/22, R200(19)/22, R200(21)/22]
b = [R200(1)/11, R200(2)/11, R200(3)/11, R200(4)/11, R200(5)/11,
     R200(6)/11, R200(7)/11, R200(8)/11, R200(9)/11, R200(10)/11]

z = R200(z0_exact)

def poch(r, n):
    """Pochhammer symbol (r)_n = gamma(r+n)/gamma(r)"""
    return gamma(r + n) / gamma(r)

def R_term(n):
    num = prod(poch(ai, n) for ai in a)
    den = prod(poch(bj, n) for bj in b)
    den *= gamma(R200(n + 1))
    return num / den

# Compute S0 = sum R(n) * z^n (the A-term)
# Compute S1 = sum R(n) * n * z^n (the B-term)
S0 = R200(0)
S1 = R200(0)
for n in range(100):
    r = R_term(n)
    zn = z**n
    term0 = r * zn
    term1 = r * n * zn
    S0 += term0
    S1 += term1
    if n > 5 and abs(term0) / abs(S0) < R200(2)**(-200):
        break

pi_val = R200(pi)
print(f"S0 = sum R(n) * z^n = {S0}")
print(f"S1 = sum R(n) * n * z^n = {S1}")
print(f"pi * S0 = {pi_val * S0}")
print(f"pi * S1 = {pi_val * S1}")

# ======================================================================
# Step 5: PSLQ search for (A,B) at exact modulus
# ======================================================================
print("\n" + "=" * 70)
print("STEP 5: PSLQ search for (A,B) at exact z0^CM")
print("=" * 70)

from sage.arith.functions import lcm

# For given (A,B), compute C as: C = pi * (A*S0 + B*S1)
# We want C^2 to be a rational number (for d=1, C should be rational*sqrt(K)/integer)

# Search over a range of A and B
print("Searching A in [1,50], B in [1,500]...")
best_error = R200(1)
best_AB = (0, 0)
best_C = None

results = []
for A in range(1, 51):
    for B in range(1, 501):
        C_pi = R200(A) * S0 + R200(B) * S1  # = C/pi
        C = pi_val * C_pi
        C2 = C * C
        
        # Try to find rational approximation of C2
        # Convert to rational with bounded denominator
        rat = C2.nearby_rational(max_denominator=50000)
        err = abs(C2 - R200(rat))
        if err < R200('1e-30'):
            C_rat = sqrt(rat)
            results.append((float(err), A, B, rat, C_rat))
            if float(err) < float(best_error):
                best_error = err
                best_AB = (A, B)
                best_C = C_rat

if results:
    results.sort(key=lambda x: x[0])
    print(f"\nTop 5 PSLQ hits at exact z0^CM:")
    for i, (err, A, B, C2_rat, C_val) in enumerate(results[:5]):
        print(f"  #{i+1}: (A,B) = ({A}, {B}), C^2 = {C2_rat}, error = {err:.2e}")
    
    print(f"\nBest candidate: (A*,B*) = {best_AB}")
    print(f"  C = {best_C}")
    print(f"  Error: {best_error:.2e}")
else:
    print("No convincing PSLQ hits found within search range.")
    print("This could mean:")
    print("  - The correct (A,B) is outside the search range")
    print("  - The parameter family does not admit a C/pi closed form")
    print("  - Higher precision is needed")

# ======================================================================
# Step 6: Compare with trial-modulus candidate (16, 250)
# ======================================================================
print("\n" + "=" * 70)
print("STEP 6: Comparison with trial-modulus candidate (16, 250)")
print("=" * 70)

A_trial = 16
B_trial = 250
C_trial = sqrt(R200(2381505)) / R200(30)

# Evaluate at exact modulus
C_pi_trial = R200(A_trial) * S0 + R200(B_trial) * S1
C_trial_at_exact = pi_val * C_pi_trial

print(f"Trial-modulus candidate (A,B) = ({A_trial}, {B_trial})")
print(f"  At exact z0^CM: C = {C_trial_at_exact}")
print(f"  Expected C (if invariant): {C_trial}")
print(f"  Relative error: {abs(C_trial_at_exact - C_trial) / abs(C_trial):.6e}")

# Check if (A,B) = (16,250) appears in our PSLQ search
found = False
for err, A, B, C2_rat, C_val in results:
    if A == 16 and B == 250:
        found = True
        print(f"\nRESULT: (A,B) = (16,250) verified at exact modulus!")
        print(f"  C = sqrt({C2_rat}) / denom")
        print(f"  Error: {err:.2e}")
        break

if not found:
    print(f"\nRESULT: (A,B) = (16,250) was NOT found in PSLQ search at exact modulus.")
    print(f"  The trial-modulus candidate does NOT survive at z0^CM.")
    if best_AB != (0, 0):
        print(f"  The correct coefficients at exact modulus appear to be (A*,B*) = {best_AB}")

print("\nDone.")
