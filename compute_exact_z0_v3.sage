"""Fix PARI polmodular parsing - vector elements are coefficients of x"""
from sage.all import *

phi = pari.polmodular(11)
print(f"phi length: {len(phi)}")

# phi[k] is coefficient of x^k in Phi_11(x, y)
# where x = j(11*tau), y = j(tau)
# Build sum_{k=0}^{12} phi[k](y) * x^k
# then substitute y = 1728

var_x = pari('x')
var_y = pari('y')

# Build polynomial in x at y = 1728
poly_x = pari(0)
for k in range(13):
    coeff = pari.subst(phi[k], var_y, 1728)
    poly_x = poly_x + coeff * var_x**k

print(f"Poly in x (degree): {pari.poldegree(poly_x)}")
print(f"Poly in x: {str(poly_x)[:200]}...")

# Now find roots
roots = pari.polroots(poly_x)
print(f"\nNumber of roots: {len(roots)}")

# Find positive real roots
pos_real_roots = []
for r in roots:
    r_real = float(pari.real(r))
    r_imag = float(pari.imag(r))
    if abs(r_imag) < 1e-6 and r_real > 0:
        pos_real_roots.append(r_real)
        print(f"  POSITIVE REAL: {r_real:.6e}")

pos_real_roots.sort()
print(f"\nPositive real roots (sorted):")
for r in pos_real_roots:
    print(f"  {r:.6e}")

# j(11i) from theta functions was ~ 1.04e30
# Let's find which root is closest to this
from sage.all import RealField
R200 = RealField(200)

# Convert j-value to lambda-value using j = 256*(1-l+l^2)^3/(l^2*(1-l)^2)
for j_val in pos_real_roots[:5]:
    jR = R200(j_val)
    # Solve for lambda: j*lambda^2*(1-lambda)^2 = 256*(1-lambda+lambda^2)^3
    # For small lambda, lambda ≈ sqrt(256/j) ≈ 16/sqrt(j)
    lam_approx = 16 / R200(j_val).sqrt()
    print(f"\nj(11i) = {jR:.6e}")
    print(f"  lambda ≈ 16/sqrt(j) = {lam_approx:.10e}")

# The first (smallest) positive root is our candidate
# Convert j-value to lambda-value properly using PARI
print("\n\nConverting each j to lambda:")
for j_val in pos_real_roots:
    jR = R200(j_val)
    # Find lambda root of j*lambda^2*(1-lambda)^2 - 256*(1-lambda+lambda^2)^3 = 0
    # Use Newton's method with lambda_0 = 16/sqrt(j)
    lam0 = 16 / jR.sqrt()
    for _ in range(10):
        l = lam0
        f = jR * l*l * (1-l)*(1-l) - 256 * (1 - l + l*l)**3
        df = jR * (2*l*(1-l)**2 - 2*l*l*(1-l)) - 256 * 3 * (1-l+l*l)**2 * (-1 + 2*l)
        lam0 = l - f/df
    if lam0 > 0 and lam0 < 1:
        print(f"  j={jR:.6e} -> lambda={float(lam0):.10e}")

print("\nDone.")
