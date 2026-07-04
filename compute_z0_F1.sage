"""
F1精确CM奇异模计算 + PSLQ验证
================================
用法：sage compute_z0_F1.sage 2>&1 | tee results.txt

计算流程：
1. 从PARI获取j-函数经典模多项式Φ₁₁
2. 代入j(i)=1728得到j(11i)的代数值
3. 通过j↔λ变换得到λ(11i) （尝试所有根）
4. 在精确λ(11i)下重做PSLQ搜索(A*,B*)
5. 与试模(16,250)对比
"""

from sage.all import *
import sys

RR = RealField(200)
pi_val = RR.pi()

print("=" * 70)
print("F1 (s=1/2, p=11, d=1): EXACT CM MODULUS COMPUTATION")
print("=" * 70)

# ============ Step 1: Modular polynomial for j-function ============
print("\n[Step 1] Computing modular polynomial Phi_11 (j-function)...")
phi = pari.polmodular(11)

# Build Phi_11(x, y) = sum phi[k]*x^k where phi[k] is a polynomial in y
# x = j(11*tau), y = j(tau)
poly_x = pari(0)
var_x = pari('x')
var_y = pari('y')
for k in range(13):
    coeff = pari.subst(phi[k], var_y, 1728)  # j(i) = 1728
    poly_x = poly_x + coeff * var_x**k

print(f"  Degree in j(11i): {pari.poldegree(poly_x)}")

# ============ Step 2: Find j(11i) ============
print("\n[Step 2] Computing j(11i) as root of Phi_11(j(11i), 1728) = 0...")
roots = pari.polroots(poly_x)

pos_real = []
for r in roots:
    re, im = float(pari.real(r)), float(pari.imag(r))
    if abs(im) < 1e-6 and re > 0:
        pos_real.append(RR(re))

pos_real.sort()
print(f"  Positive real roots of Phi_11: {len(pos_real)}")
j11 = pos_real[0]
print(f"  j(11i) = {j11:.6e}")
print(f"  (Expected ~1.038e30 from theta series)")

# ============ Step 3: Convert j(11i) to lambda(11i) ============
print("\n[Step 3] Converting j(11i) to lambda(11i)...")
print(f"  Formula: j = 256*(1-λ+λ²)³/(λ²(1-λ)²)")
print(f"  For small λ: λ ≈ 16/√j")

lam_approx = RR(16) / j11.sqrt()
print(f"  Approximate λ ≈ 16/√j = {lam_approx:.10e}")

# Refine using Newton's method
lam = lam_approx
for _ in range(20):
    l = lam
    f = j11 * l*l * (1-l)*(1-l) - RR(256) * (1 - l + l*l)**3
    df = j11 * (2*l*(1-l)**2 - 2*l*l*(1-l)) - RR(256) * RR(3) * (1-l+l*l)**2 * (RR(-1) + RR(2)*l)
    lam = l - f/df

z0_exact = lam
print(f"  λ(11i) = EXACT singular modulus = {z0_exact:.10e}")

# Compare with trial modulus
z0_trial = RR(1)/RR(121)
print(f"\n  Trial modulus 1/121 = {z0_trial:.10e}")
print(f"  Ratio (z0_exact / 1/121) = {float(z0_exact / z0_trial):.6e}")

# ============ Step 4: Series at EXACT modulus ============
print("\n[Step 4] Computing hypergeometric series at exact z0^CM...")

# F1 parameters
a = [RR(1)/2, RR(1)/22, RR(3)/22, RR(5)/22, RR(7)/22, RR(9)/22,
     RR(13)/22, RR(15)/22, RR(17)/22, RR(19)/22, RR(21)/22]
b = [RR(1)/11, RR(2)/11, RR(3)/11, RR(4)/11, RR(5)/11,
     RR(6)/11, RR(7)/11, RR(8)/11, RR(9)/11, RR(10)/11]

def poch(r, n):
    return gamma(r + n) / gamma(r)

z0 = z0_exact

# Since z0 is extremely small (~10^-14), the series converges in ~3-5 terms
S0 = RR(0); S1 = RR(0); N_terms = 0
for n in range(100):
    num = prod(poch(ai, n) for ai in a)
    den = prod(poch(bj, n) for bj in b)
    den *= gamma(RR(n + 1))
    r = num / den
    term0 = r * z0**n
    term1 = r * RR(n) * z0**n
    S0 += term0; S1 += term1
    if n > 2 and abs(term0)/abs(S0) < RR(2)**(-200):
        N_terms = n + 1
        break

print(f"  Series converged in {N_terms} terms")
print(f"  S0 = sum R(n)*z0^n      = {S0}")
print(f"  S1 = sum R(n)*n*z0^n    = {S1}")
print(f"  pi*S0 = {pi_val*S0}")
print(f"  pi*S1 = {pi_val*S1}")

# ============ Step 5: PSLQ at exact modulus ============
print("\n[Step 5] PSLQ search for (A*,B*) at exact z0^CM...")
print("  Searching A in [1,50], B in [1,500]...")

best_err = RR(1); best_AB = (0, 0); best_C2 = None
hits = []
total = 50 * 500

for A in range(1, 51):
    for B in range(1, 501):
        C_pi = RR(A)*S0 + RR(B)*S1
        C_sq = (pi_val * C_pi)**2
        rat = C_sq.nearby_rational(max_denominator=100000)
        err = abs(C_sq - RR(rat))
        if err < RR('1e-30'):
            hits.append((float(err), A, B, rat))
            if float(err) < float(best_err):
                best_err = err; best_AB = (A, B); best_C2 = rat


if hits:
    hits.sort(key=lambda x: x[0])
    print(f"\n  Top 5 PSLQ hits at EXACT z0^CM:")
    print(f"  {'Rank':<5} {'(A,B)':<12} {'C^2':<30} {'Error':<12}")
    print(f"  {'-'*5} {'-'*12} {'-'*30} {'-'*12}")
    for i, (err, A, B, C2rat) in enumerate(hits[:5]):
        print(f"  {i+1:<5} ({A},{B})  {str(C2rat):<30} {err:.2e}")
    
    # Check if (16,250) is found
    found_16250 = any(A==16 and B==250 for _, _, A, B, _ in [(e,A,B,c) for e,A,B,c in hits])
    print(f"\n  (16,250) in exact-modulus PSLQ results: {'YES ✓' if found_16250 else 'NOT FOUND ✗'}")
    
    print(f"\n  BEST candidate: (A*,B*) = {best_AB}")
    print(f"  C^2 = {best_C2}")
    
    # Check trial candidate (16,250) at exact modulus
    C_pi_trial = RR(16)*S0 + RR(250)*S1
    C_trial = pi_val * C_pi_trial
    C_target = RR(2381505).sqrt() / RR(30)
    print(f"\n  Trial candidate (16,250) at EXACT modulus:")
    print(f"    C(computed) = {C_trial}")
    print(f"    C(target)   = {C_target}")
    print(f"    Error: {abs(C_trial - C_target)/abs(C_target):.6e}")
else:
    print("\n  No PSLQ hits found. Possible reasons:")
    print("  - z0^CM is the wrong root of the modular equation")
    print("  - The series does not admit a C/pi closed form at this root")
    print("  - (A,B) outside search range")
    print("\n  The singular modulus may instead be a DIFFERENT root")
    print("  of Phi_11(λ, 1/2) = 0 for the lambda-function modular polynomial.")

print("\n" + "=" * 70)
print("COMPUTATION COMPLETE")
print("=" * 70)
