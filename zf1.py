"""F1 conjecture: SageMath 200-bit verification"""
from sage.all import *

RF = RealField(200)
O = RF(1)

# F1: s=1/2, p=11 — all parameters as RF floats
a = [O/2, O/22, RF(3)/22, RF(5)/22, RF(7)/22, RF(9)/22,
     RF(13)/22, RF(15)/22, RF(17)/22, RF(19)/22, RF(21)/22]
b = [O/11, RF(2)/11, RF(3)/11, RF(4)/11, RF(5)/11,
     RF(6)/11, RF(7)/11, RF(8)/11, RF(9)/11, RF(10)/11]

def Rterm(n):
    num = prod(gamma(ai + n) / gamma(ai) for ai in a)
    den = prod(gamma(bj + n) / gamma(bj) for bj in b)
    den *= gamma(RF(n + 1))
    return num / den

A = RF(16); B = RF(250)
z = O / RF(121)

S = RF(0)
for n in range(100):
    t = Rterm(n) * (A + B * n) * z**n
    S += t
    if n > 5 and abs(t) / abs(S) < RF(10)**(-80):
        break

target = RF(2381505).sqrt() / (RF(30) * RF.pi())
err = abs(S - target) / abs(target)

print(f"SageMath 10.9")
print(f"Sum R(n)(16+250n)/121^n = sqrt(2381505)/(30*pi)")
print(f"  Sum:    {S}")
print(f"  Target: {target}")
print(f"  Error:  {err}")
print(f"  Status: {'PASSES' if err < 1e-30 else 'FAILS'}")
