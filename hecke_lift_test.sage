"""Compute Hecke lift T_p(lambda)(i) = (1/p)*sum_{ad=p, 0<=b<d} lambda((a*tau+b)/d)
This is the correct singular modulus for the descent operator.
"""
from sage.all import *

def lam_via_theta(tau, dps=500):
    """Compute lambda(tau) using theta functions at high precision"""
    RR = RealField(dps)
    old_dps = RR.prec()
    
    # q = exp(pi*i*tau)  
    q = exp(pi * ComplexField(dps).I() * ComplexField(dps)(tau))
    
    # theta_2 and theta_3
    t2 = ComplexField(dps).zero()
    t3 = ComplexField(dps).zero()
    for n in range(200):
        term2 = q**(n*(n+1))
        t2 += term2
        term3 = q**(n*n) if n > 0 else ComplexField(dps).one()
        t3 += 2*term3 if n > 0 else ComplexField(dps).one()
        if abs(term2) < RR(2)**(-dps) and abs(term3) < RR(2)**(-dps) if n > 0 else False:
            break
    
    t2 = 2 * q**(RR(1)/4) * t2
    lam = (t2**4) / (t3**4)
    return lam.real()

# Test: Hecke lift for p=7, should equal 1/2401 if this is the right approach
print("Testing Hecke lift T_7(lambda)(i):")
RR200 = RealField(200)
tau0 = RR200(0) + RR200(1)*ComplexField(200).I()  # tau = i

# T_p f(tau) = (1/p) * sum_{ad=p, 0<=b<d} f((a*tau+b)/d)
# For prime p: a in {1, p}, d in {p, 1}
# a=1, d=p: terms (tau + b)/p for b=0,...,p-1
# a=p, d=1: term p*tau

# For p=7, compute T_7(lambda)(i) with low precision first
p = 7
dps = 100
CC = ComplexField(dps)
I = CC.I()
tau = I

# Hecke lift terms
terms = []
# a=p, d=1: f(p*tau)
terms.append((p*tau).real())  # This is 7i, lambda(7i) is real
# a=1, d=p: f((tau+b)/p) for b=0,...,p-1
for b in range(p):
    terms.append((tau + CC(b)) / p)

print(f"Number of terms in Hecke sum: {len(terms)}")
print(f"Terms (argument):")
for i, t in enumerate(terms):
    print(f"  term[{i}]: tau = {t}")
"""
