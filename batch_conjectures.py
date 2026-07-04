"""
Batch generate numerical conjectures for F2, F5, F9, F13.
Same method as F1: trial z0=1/p^2, integer search for (A,B), rational C^2.
"""
from mpmath import mp, mpf, sqrt, pi, gamma
from fractions import Fraction

mp.dps = 300

ATOMIC = {
    "s=1/2": {"a": [Fraction(1,2), Fraction(1,2)], "b": [Fraction(1,1)], "N": 2},
    "s=1/3": {"a": [Fraction(1,3), Fraction(2,3)], "b": [Fraction(1,1)], "N": 3},
    "s=1/4": {"a": [Fraction(1,4), Fraction(3,4)], "b": [Fraction(1,1)], "N": 4},
    "s=1/6": {"a": [Fraction(1,6), Fraction(5,6)], "b": [Fraction(1,1)], "N": 6},
}

def translation_Sp(seed_a, seed_b, p):
    def split(r):
        return [Fraction(r.numerator + j*r.denominator, p*r.denominator) for j in range(p)]
    def norm(params):
        seen, result = set(), []
        for r in params:
            red = Fraction(r.numerator % r.denominator, r.denominator)
            if red.numerator == 0 or red == Fraction(1,1): continue
            key = (red.numerator, red.denominator)
            if key not in seen: seen.add(key); result.append(red)
        return sorted(result, key=lambda x: (x.denominator, x.numerator))
    na = []; [na.extend(split(ai)) for ai in seed_a]; na = norm(na)
    nb = []; [nb.extend(split(bj)) for bj in seed_b]
    for j in range(1, p): nb.append(Fraction(j, p))
    nb = norm(nb)
    return na, nb

def poch(r, n):
    x = mpf(r.numerator) / r.denominator
    return gamma(x + n) / gamma(x)

def series_term(a, b, z0, A, B, max_n=200):
    """Sum R(n)(A+Bn)z0^n"""
    S = mpf(0)
    for n in range(max_n):
        num = mpf(1)
        for ai in a: num *= poch(ai, n)
        den = mpf(1)
        for bj in b: den *= poch(bj, n)
        for _ in range(len(a) - len(b)): den *= gamma(n+1)
        rn = num / den * z0**n
        S += rn * (A + B*n)
        if n > 5 and abs(rn * (A + B*n)) / max(abs(S), mpf('1e-100')) < mpf('1e-60'):
            break
    return S, n

def find_ABC(a, b, p, family_name):
    """Find (A,B) and C for given hypergeometric parameters at trial z0=1/p^2"""
    z0 = mpf(1) / mpf(p*p)
    print(f"\n{'='*65}")
    print(f"FAMILY {family_name}: s (see atomic), p={p}")
    print(f"  |a|={len(a)}, |b|={len(b)}, k={len(a)-len(b)}")
    print(f"  z0 = 1/{p**2} = {mp.nstr(z0, 12)}")

    # Compute S0, S1
    S0, S1 = mpf(0), mpf(0)
    for n in range(200):
        num = mpf(1)
        for ai in a: num *= poch(ai, n)
        den = mpf(1)
        for bj in b: den *= poch(bj, n)
        for _ in range(len(a) - len(b)): den *= gamma(n+1)
        rn = num / den * z0**n
        S0 += rn; S1 += rn * n
        if n > 5 and abs(rn) * (1 + 200*n) / max(abs(S0), 1e-100) < mpf('1e-80'):
            break

    # Search (A,B): (A*S0 + B*S1)*pi should be algebraic
    best = None
    for A_try in range(1, 30):
        for B_try in range(1, 301):
            C_val = pi * (A_try * S0 + B_try * S1)
            C_sq = float(C_val * C_val)
            # Is C^2 rational? Check denominators up to 500
            for den in range(1, 501):
                num = round(C_sq * den)
                err = abs(C_sq - num/den)
                if err < 0.001 and num > 0:
                    quality = err * den  # weighted by complexity
                    if best is None or err < best['err']:
                        best = {'A': A_try, 'B': B_try, 'C_sq_num': num, 'C_sq_den': den,
                                'C_val': float(C_val), 'err': err}

    if best:
        N, D = best['C_sq_num'], best['C_sq_den']
        # Try to express C^2 = N*K / D^2 (standard form C = sqrt(NK)/D)
        best_form = None
        for nn in range(1, 200):
            for dd in range(1, 200):
                target = N / D
                KK = round(target * dd * dd / (nn * nn))
                for kk in [KK-1, KK, KK+1]:
                    if kk <= 0: continue
                    predicted = nn * nn * kk / (dd * dd)
                    err2 = abs(target - predicted)
                    if err2 < 0.01 and (best_form is None or err2 < best_form['err']):
                        best_form = {'N': nn, 'D': dd, 'K': kk, 'err': err2}

        print(f"  Best (A,B): ({best['A']}, {best['B']})")
        print(f"  C ≈ {best['C_val']:.10f}")
        print(f"  C^2 ≈ {best['C_sq_num']}/{best['C_sq_den']} (err={best['err']:.6f})")

        if best_form:
            NF = best_form
            C_exact = mpf(NF['N']) * sqrt(mpf(NF['K'])) / mpf(NF['D'])
            print(f"  Standard form: C = {NF['N']}*sqrt({NF['K']})/{NF['D']} = {mp.nstr(C_exact, 15)}")
            # Verify
            true_S, terms = series_term(a, b, z0, best['A'], best['B'])
            target = C_exact / pi
            rel_err = abs(true_S - target) / abs(target)
            print(f"  Sum[{best['A']}+{best['B']}n] = {mp.nstr(true_S, 15)}")
            print(f"  Target C/pi = {mp.nstr(target, 15)}")
            print(f"  Rel error = {mp.nstr(rel_err, 8)}")
            print(f"  [{'PASS' if rel_err < 1e-10 else 'CLOSE'}]")

            # Output LaTeX conjecture
            print(f"\n  === LATEX ===")
            print(f"  \\\\begin{{conjecture}}[{family_name} family]")
            print(f"  At trial $z_0=1/{p**2}$, the hypergeometric series")
            print(f"  $\\\\sum R(n)({best['A']}+{best['B']}n)/{p**2}n \\\\approx"
                  f" {NF['N']}\\\\sqrt{{{NF['K']}}}/({NF['D']}\\\\pi)$")
            print(f"  to $\\\\sim{mp.nstr(rel_err, 6)}$ relative error at 300-digit precision.")
            print(f"  \\\\end{{conjecture}}")
    else:
        print(f"  No clean (A,B) found at this z0. Try different p or z0 pattern.")

# Run for all 4 families
families = [
    ("F2", "s=1/2", 13),
    ("F5", "s=1/3", 11),
    ("F9", "s=1/4", 11),
    ("F13", "s=1/6", 11),
]

for fam, seed, p in families:
    sd = ATOMIC[seed]
    a, b = translation_Sp(sd["a"], sd["b"], p)
    find_ABC(a, b, p, fam)
