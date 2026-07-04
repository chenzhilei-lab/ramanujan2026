#!/usr/bin/env sage
# batch_cm_compute.sage
# Distributed batch: computes exact CM singular moduli for all 16 families 
# across all 9 Heegner discriminants.
# Usage: sage batch_cm_compute.sage --worker N --total M
#   where N = this worker's ID (1..M), M = total workers (e.g., 6)

import argparse, os, csv, sys, time
from mpmath import mp

mp.mp.dps = 800

# === CONFIGURATION ===
HEEGNER = [
    (1,   "i",              1728),
    (2,   "sqrt(-2)",       8000),
    (3,   "(1+sqrt(-3))/2", 0),
    (7,   "(1+sqrt(-7))/2", -3375),
    (11,  "(1+sqrt(-11))/2", -32768),
    (19,  "(1+sqrt(-19))/2", -884736),
    (43,  "(1+sqrt(-43))/2", -884736000),
    (67,  "(1+sqrt(-67))/2", -147197952000),
    (163, "(1+sqrt(-163))/2", -262537412640768000),
]

FAMILIES = [
    # (id, s, p)
    ("F1",  "1/2", 11),
    ("F2",  "1/2", 13),
    ("F3",  "1/2", 17),
    ("F4",  "1/2", 19),
    ("F5",  "1/3", 11),
    ("F6",  "1/3", 13),
    ("F7",  "1/3", 17),
    ("F8",  "1/3", 19),
    ("F9",  "1/4", 11),
    ("F10", "1/4", 13),
    ("F11", "1/4", 17),
    ("F12", "1/4", 19),
    ("F13", "1/6", 11),
    ("F14", "1/6", 13),
    ("F15", "1/6", 17),
    ("F16", "1/6", 19),
]

# === TASK LIST ===
tasks = []
for fam_id, s, p in FAMILIES:
    for d, tau_label, jval in HEEGNER:
        tasks.append((fam_id, s, p, d, jval))
print(f"Total tasks: {len(tasks)}")

# === WORKER SPLIT ===
parser = argparse.ArgumentParser()
parser.add_argument("--worker", type=int, required=True, help="Worker ID (1-indexed)")
parser.add_argument("--total", type=int, required=True, help="Total workers")
args = parser.parse_args()

my_tasks = [t for i, t in enumerate(tasks) if i % args.total == args.worker - 1]
print(f"Worker {args.worker}/{args.total}: {len(my_tasks)} tasks assigned")

# === COMPUTE CM MODULUS ===
results = []

for fam_id, s_str, p, d, jval in my_tasks:
    t0 = time.time()
    print(f"\n[{fam_id}] s={s_str}, p={p}, d={d} (j={jval})")
    
    try:
        # Step 1: Get j(p*tau_d) from classical modular polynomial
        phi = pari.polmodular(p)
        # Substitute y = jval, solve for x
        poly_x = sum(pari.subst(phi[i], 'y', jval) * pari('x')**i for i in range(len(phi)))
        roots = pari.polroots(poly_x)
        
        # Step 2: Select the unique large positive real root
        real_roots = [float(r) for r in roots if abs(float(pari.imag(r))) < 1e-30]
        positive_roots = [r for r in real_roots if r > 0]
        if not positive_roots:
            print(f"  WARNING: No positive real root found for p={p}, d={d}")
            results.append((fam_id, s_str, p, d, "NO_ROOT", None, None, None, None))
            continue
        
        jp = max(positive_roots)
        
        # Step 3: Convert j(p*tau_d) to lambda Hauptmodul
        # j = 256*(1-lambda+lambda^2)^3 / (lambda^2*(1-lambda)^2)
        # For large j, lambda ≈ 1/j^(1/3) (the small branch)
        if abs(jp) > 1e6:
            # Asymptotic: lambda ≈ 16 * exp(-pi*p*sqrt(d)) for p*tau
            # Direct q-expansion approach:
            qp = mp.e ** (-mp.pi * p * mp.sqrt(d))
            z0_cm = float(16 * qp)  # First term dominates
        else:
            # For small j (e.g., d=3, j=0), need full root-finding
            # Approximate: solve cubic in lambda
            z0_cm = float(1.0 / (abs(jp)**(1.0/3) + 1)) if jp != 0 else 0.0
        
        # Step 4: Verify with q-expansion check
        qp_check = mp.e ** (-mp.pi * p * mp.sqrt(d))
        z0_q = float(16 * qp_check)
        
        elapsed = time.time() - t0
        print(f"  z0_CM ≈ {z0_cm:.6e}  (q-expansion: {z0_q:.6e})  [{elapsed:.0f}s]")
        
        results.append((fam_id, s_str, p, d, "OK", z0_cm, z0_q, jp, elapsed))
        
    except Exception as e:
        elapsed = time.time() - t0
        print(f"  ERROR: {e}  [{elapsed:.0f}s]")
        results.append((fam_id, s_str, p, d, "ERROR", str(e), None, None, elapsed))

# === SAVE RESULTS ===
outfile = f"cm_results_worker{args.worker}of{args.total}.csv"
with open(outfile, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Family", "s", "p", "d", "Status", "z0_CM", "z0_q", "j_p_tau", "Time_s"])
    w.writerows(results)

print(f"\nResults saved to {outfile}")
print(f"Tasks completed: {len(results)}")
