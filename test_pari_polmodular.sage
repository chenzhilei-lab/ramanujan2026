"""Fix variable conventions for PARI polmodular"""
from sage.all import *

phi = pari.polmodular(11)
print("phi type:", type(phi))
print("phi length:", len(phi))

f = phi[0]
print("\nVariables in polynomial:", pari.variables(f))
print("\nPolynomial (first 200 chars):", str(f)[:200])

# In PARI polmodular(N, {X}, {Y}):
# The polynomial Phi_N(X, Y) where X = j(N*tau), Y = j(tau)
# polmodular(N) returns F where F(j(N*tau), j(tau)) = 0
# By default: X is 'x', Y is 'y' (or vice versa)

# Substitute Y = j(i) = 1728
g = pari.subst(f, 'y', 1728)
print("\nAfter subst(y, 1728):", str(g)[:200], "...")
print("Degree in x:", pari.poldegree(g, 'x'))

# Find all real roots
R = PolynomialRing(QQ, 'x')
x = R.gen()
G = R(pari.Vecrev(pari.polcoeff(g, 0)) if pari.type(g) == 't_POL' else g)
# Actually, g is a polynomial in x
# Let me convert properly
G = R(str(g).replace('y', ''))  # This is hacky, let me do it properly

# Actually, after substituting y=1728, g is a polynomial in x
# with huge coefficients, but in PARI representation
# Let me use a different approach
RR = RealField(200)
# Evaluate g at various points to bracket roots
for log_val in [20, 25, 28, 30, 32, 35]:
    val = RR(10)**log_val
    result = pari.subst(g, 'x', val)
    # Convert to RR
    try:
        r = RR(str(result))
        print(f"  g({val:.0e}) = {float(r):.2e}")
    except:
        pass

# Try to find the root using Newton's method starting from 1e30
print("\nTrying PARI's polroots on g...")
roots = pari.polroots(g)
print(f"Number of roots: {len(roots)}")
for i in range(min(5, len(roots))):
    r = roots[i]
    r_real = pari.real(r)
    r_imag = pari.imag(r)
    print(f"  root[{i}]: {float(r_real):.6e} + {float(r_imag):.6e}i")

# Find the positive real root closest to 1e30
target = 1.03819703567651e30
best_root = None
best_dist = float('inf')
for r in roots:
    r_real = float(pari.real(r))
    r_imag = float(pari.imag(r))
    if abs(r_imag) < 1e-10 and r_real > 0:
        dist = abs(r_real - target)
        if dist < best_dist:
            best_dist = dist
            best_root = r_real

print(f"\nBest root for j(11i) = {best_root:.6e}")
print(f"Expected ~ 1.038e30, ratio = {best_root/target:.6f}")
