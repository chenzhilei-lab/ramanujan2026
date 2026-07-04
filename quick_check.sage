"""Quick check: does (16,250) work at trial modulus with 500-bit precision?"""
from sage.all import *
RR = RealField(500)
pi_val = RR.pi()

a = [RR(1)/2, RR(1)/22, RR(3)/22, RR(5)/22, RR(7)/22, RR(9)/22,
     RR(13)/22, RR(15)/22, RR(17)/22, RR(19)/22, RR(21)/22]
b = [RR(1)/11, RR(2)/11, RR(3)/11, RR(4)/11, RR(5)/11,
     RR(6)/11, RR(7)/11, RR(8)/11, RR(9)/11, RR(10)/11]

def poch(r, n):
    return gamma(r + n) / gamma(r)

z0 = RR(1)/121
S0 = RR(0); S1 = RR(0)
for n in range(100):
    num = prod(poch(ai, n) for ai in a)
    den = prod(poch(bj, n) for bj in b)
    den *= gamma(RR(n + 1))
    r = num / den
    S0 += r * z0**n
    S1 += r * RR(n) * z0**n
    if n > 2 and abs(r*z0**n)/abs(S0) < RR(2)**(-500): break

pS0 = pi_val * S0
pS1 = pi_val * S1
print(f"pi*S0 = {pS0}")
print(f"pi*S1 = {pS1}")

# Check (16,250):
C = RR(16)*pS0 + RR(250)*pS1
C2 = C**2
print(f"\n(16,250): C/pi = {C}")
print(f"C = {pi_val*C:.20e}")
print(f"C^2 = {C2}")
print(f"Expected C^2 = 2381505/900 = {RR(2381505)/900}")
print(f"Error = {abs(C2 - RR(2381505)/900):.2e}")

# Check nearby rational
rat = C2.nearby_rational(max_denominator=100000)
print(f"Nearby rational = {rat}")
print(f"Error = {abs(C2 - RR(rat)):.2e}")
