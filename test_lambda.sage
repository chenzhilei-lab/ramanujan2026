"""Compute lambda(pi) using PARI's modular forms"""
from sage.all import pari

# PARI's elllambda function computes the modular lambda function
# elllambda(tau, n) returns the values lambda and 1-lambda at n*tau
# Actually, in PARI, the function is different...

# Let me try using PARI's theta functions directly
print("Computing lambda(pi) for p=2,3,5,7,11:")
for p in [2, 3, 5, 7, 11]:
    tau = pari(p) * pari.I()
    # Compute q = exp(I*pi*tau) = exp(-p*pi)
    q = pari.exp(pari.I() * pari.pi() * tau)
    
    # theta_2 = 2*q^(1/4) * sum_{n>=0} q^(n*(n+1))
    t2 = 2 * q**0.25
    s = pari.zero()
    for n in range(100):
        term = q**(n*(n+1))
        s += term
        if abs(term) < pari('10^-100'):
            break
    t2 *= s
    
    # theta_3 = 1 + 2*sum q^(n^2)
    t3 = pari(1)
    for n in range(1, 100):
        term = q**(n*n)
        t3 += 2 * term
        if abs(term) < pari('10^-100'):
            break
    
    lam = (t2**4 / t3**4).real()
    print(f"  p={p}: lambda({p}i) = {lam}")
    
    # Also compute the j-invariant at p*i
    j = pari.ellj(tau)
    print(f"        j({p}i) = {j}")

# Verify: j(i) should be 1728
print(f"\nj(i) = {pari.ellj(pari.I())}")

# For p=2: j(2i) should be related to (sqrt(2)-1)^4
print(f"j(2i) = {pari.ellj(2*pari.I())}")
print(f"(sqrt(2)-1)^4 = {(pari.sqrt(pari(2))-1)**4}")
