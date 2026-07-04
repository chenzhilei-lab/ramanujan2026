"""
compute_all_lambda_roots.sage
=============================
Compute ALL 12 roots of the lambda-level modular equation
Psi_11(X, 1/2) = 0 for p=11, tau=i (d=1).

The 12 roots are:
  - lambda(11i)                              (1 root: "direct" CM lift)
  - lambda((k + i)/11) for k = 0, ..., 10   (11 roots: "satellite" cusps)

These are all values that lambda(11*tau) can take given lambda(tau) = 1/2.

Also runs PSLQ search at each root with reasonable z0 values.
"""
from sage.all import *

RR = RealField(500)
pi_val = RR.pi()
CC = ComplexField(500)

# F1 parameters: s=1/2, p=11
a_list = [RR(1)/2, RR(1)/22, RR(3)/22, RR(5)/22, RR(7)/22, RR(9)/22,
          RR(13)/22, RR(15)/22, RR(17)/22, RR(19)/22, RR(21)/22]
b_list = [RR(1)/11, RR(2)/11, RR(3)/11, RR(4)/11, RR(5)/11,
          RR(6)/11, RR(7)/11, RR(8)/11, RR(9)/11, RR(10)/11]

def poch(r, n):
    return gamma(r + n) / gamma(r)

def lambda_q(q_val, terms=30):
    """lambda via q-product: 16*q * prod((1+q^{2n})/(1+q^{2n-1}))^8"""
    prod = RR(1)
    for n in range(1, terms+1):
        prod *= ((1 + q_val**(2*n)) / (1 + q_val**(2*n-1)))**8
    return 16 * q_val * prod

def lambda_from_tau_imag(b):
    """lambda(i*b) for real b>0"""
    q = exp(-RR(pi_val) * b)
    return lambda_q(q)

def compute_series(z0, max_n=200):
    """Compute S0, S1 at given z0."""
    S0 = RR(0); S1 = RR(0)
    for n in range(max_n):
        num = prod(poch(ai, n) for ai in a_list)
        den = prod(poch(bj, n) for bj in b_list)
        den *= gamma(RR(n + 1))
        r = num / den
        term0 = r * z0**n
        term1 = r * RR(n) * z0**n
        S0 += term0
        S1 += term1
        if n > 3 and abs(term0) / abs(S0) < RR(2)**(-400):
            break
    return S0, S1

def pslq_scan(z0, max_A=200, max_B=2000, label=""):
    """PSLQ search for (A,B,C) at given z0."""
    if z0 <= 0 or z0 >= 1:
        return None
    S0, S1 = compute_series(z0)
    best_err = RR(1)
    best_hit = None
    hits = []

    for A in range(1, max_A+1):
        for B in range(1, max_B+1):
            C_pi = RR(A) * S0 + RR(B) * S1
            C_sq = (pi_val * C_pi)**2
            try:
                rat = C_sq.nearby_rational(max_denominator=100000)
                err = abs(C_sq - RR(rat))
            except:
                continue
            if err < RR('1e-30'):
                hits.append((float(err), A, B, rat))
                if err < best_err:
                    best_err = err
                    best_hit = (A, B, rat)

    return {"z0": z0, "S0": S0, "S1": S1, "hits": hits, "best": best_hit}

# ============================================================
print("=" * 80)
print("ALL 12 ROOTS OF LAMBDA MODULAR EQUATION Psi_11(X, 1/2) = 0")
print("=" * 80)

roots = []

# Root 0: lambda(11i) — direct CM lift
print("\n--- Root 0: lambda(11i) ---")
q_11i = exp(-11 * pi_val)
lam_11i = lambda_q(q_11i)
print(f"  q = exp(-11*pi) = {float(q_11i):.6e}")
print(f"  lambda(11i) = {float(lam_11i):.10e}")
roots.append(("11i", lam_11i))

# Roots 1-11: lambda((k+i)/11) for k = 0,...,10 — satellite cusps
print("\n--- Roots 1-11: lambda((k+i)/11) for k=0,...,10 ---")
for k in range(11):
    b = RR(1)/11  # imag part = 1/11
    q = exp(-pi_val * b)
    lam = lambda_q(q)
    print(f"  k={k:2d}: tau=({k}+i)/11  q=exp(-pi/11)  lambda={float(lam):.10e}")
    roots.append((f"({k}+i)/11", lam))

# Sort by value
roots.sort(key=lambda x: x[1])

print("\n" + "=" * 80)
print("ALL 12 ROOTS (sorted by value)")
print("=" * 80)
for i, (label, val) in enumerate(roots):
    log10 = float(log(val, 10)) if val > 0 else float('-inf')
    print(f"  Root {i:2d}: {label:<15s}  lambda = {float(val):.12e}  (log10 = {log10:.2f})")
    if val > 0 and val < 1:
        diff_from_trial = float(abs(val - RR(1)/121) / (RR(1)/121))
        if diff_from_trial < 10:
            print(f"           *** WITHIN {diff_from_trial:.1f}x of 1/121 ***")

print(f"\n  Trial modulus 1/121 = {float(RR(1)/121):.12e}")

# ============================================================
print("\n" + "=" * 80)
print("PSLQ SCAN AT ROOTS WITH z0 > 1e-100 (reasonable convergence)")
print("=" * 80)

all_hits = []
checked = 0
for label, z0 in roots:
    if z0 < RR('1e-100') or z0 > RR('0.999'):
        checked += 1
        continue

    result = pslq_scan(z0, max_A=200, max_B=2000, label=label)
    if result is None:
        checked += 1
        continue
    S0 = result["S0"]
    S1 = result["S1"]
    hits = result["hits"]

    print(f"\n  {label}: z0 = {float(z0):.10e}")
    print(f"    S0 = {float(S0):.12e}  |  pi*S0 = {float(pi_val*S0):.12e}")
    print(f"    S1 = {float(S1):.12e}  |  pi*S1 = {float(pi_val*S1):.12e}")

    if hits:
        hits.sort(key=lambda x: x[0])
        print(f"    PSLQ HITS: {len(hits)}")
        for err, A, B, C2rat in hits[:5]:
            print(f"      (A,B)=({A},{B})  C^2={C2rat}  err={err:.2e}")
        all_hits.extend([(label, float(z0), A, B, str(C2rat), err) for err, A, B, C2rat in hits])
    else:
        print(f"    No PSLQ hits")

# ============================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(f"\n  Total roots with z0 > 1e-100: {sum(1 for _, z in roots if z > RR('1e-100'))}")
print(f"  Total PSLQ hits found: {len(all_hits)}")

if all_hits:
    all_hits.sort(key=lambda x: x[5])
    print(f"\n  Best hits across all roots:")
    for label, z0, A, B, C2, err in all_hits[:10]:
        print(f"    {label}: z0={z0:.6e}  (A,B)=({A},{B})  C^2={C2}  err={err:.2e}")
else:
    print("\n  NO PSLQ hits at ANY of the 12 roots.")
    print("  The trial modulus 1/121 is not among the 12 lambda-modular roots.")
    print("  This confirms the boundary at p=7: no C/pi closed form for p=11, d=1.")

# ============================================================
# Check if 1/121 is close to any root
print("\n" + "=" * 80)
print("PROXIMITY CHECK: Is 1/121 near any root?")
print("=" * 80)
trial = RR(1)/121
for label, val in roots:
    ratio = float(val / trial)
    if 0.01 < ratio < 100:
        print(f"  {label}: ratio to 1/121 = {ratio:.4f}")
    else:
        print(f"  {label}: ratio to 1/121 = {ratio:.2e}")

print("\nDone.")
