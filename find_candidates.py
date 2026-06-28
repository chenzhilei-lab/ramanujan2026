# Supplementary Material for paper1_v8_12
#
# This Python script implements the descent operator S_p as defined in
# "A Combinatorial Descent Operator for Systematic Enumeration of 
#  Admissible Hypergeometric Parameters in Ramanujan-Type 1/pi Series"
#
# Usage: python3 find_candidates.py > full_catalog.csv

from fractions import Fraction
from math import gcd

def frac(r):
    """Reduce a rational to its representative in (0,1] modulo Z."""
    if r.denominator == 0:
        return None
    num = r.numerator % r.denominator
    if num == 0:
        return None  # discard integers
    return Fraction(num, r.denominator)

def norm(param_set):
    """Apply norm(): reduce mod Z, discard 0/1, deduplicate, sort."""
    reduced = []
    seen = set()
    for r in param_set:
        f = frac(r)
        if f is not None and f not in seen:
            reduced.append(f)
            seen.add(f)
    return sorted(reduced, key=lambda x: (x.denominator, x.numerator))

def descent_operator(seed_a, seed_b, p):
    """Apply S_p to an atomic seed (a, b) at prime p."""
    # Split numerator parameters
    a_raw = []
    for ai in seed_a:
        for j in range(p):
            a_raw.append(Fraction(ai.numerator + j * ai.denominator, ai.denominator * p))
    
    # Split denominator parameters + factorial splitting
    b_raw = []
    for bi in seed_b:
        for j in range(p):
            b_raw.append(Fraction(bi.numerator + j * bi.denominator, bi.denominator * p))
    # Factorial splitting: (1)_{pn} adds {1/p, 2/p, ..., (p-1)/p}
    for j in range(1, p):
        b_raw.append(Fraction(j, p))
    
    return norm(a_raw), norm(b_raw)

def admissible_primes(N_s, pmax=31):
    """Return all primes p <= pmax with gcd(p, N_s) = 1."""
    primes = [p for p in range(2, pmax+1) 
              if all(p % d != 0 for d in range(2, int(p**0.5)+1))]
    return [p for p in primes if gcd(p, N_s) == 1]

# The four atomic seeds
seeds = {
    '1/2': (Fraction(1,2), Fraction(1,2)),  # a = (1/2, 1/2), b = (1)
    '1/3': (Fraction(1,3), Fraction(2,3)),  # a = (1/3, 2/3), b = (1)
    '1/4': (Fraction(1,4), Fraction(3,4)),  # a = (1/4, 3/4), b = (1)
    '1/6': (Fraction(1,6), Fraction(5,6)),  # a = (1/6, 5/6), b = (1)
}
N_s_values = {'1/2': 2, '1/3': 3, '1/4': 2, '1/6': 6}

# Heegner discriminants (class number 1)
heegner_d = [1, 2, 3, 7, 11, 19, 43, 67, 163]

print("s,p,d,|a'|,|b'|,k,a',b'")
total = 0
for s, (a1, a2) in seeds.items():
    N_s = N_s_values[s]
    primes = admissible_primes(N_s)
    for p in primes:
        a_prime, b_prime = descent_operator([a1, a2], [Fraction(1,1)], p)
        k = len(a_prime) - len(b_prime)
        for d in heegner_d:
            a_str = ";".join(str(x) for x in a_prime)
            b_str = ";".join(str(x) for x in b_prime)
            print(f"{s},{p},{d},{len(a_prime)},{len(b_prime)},{k},{a_str},{b_str}")
            total += 1

# Verify: should print 351
import sys
print(f"# Total: {total} structurally admissible triples", file=sys.stderr)
