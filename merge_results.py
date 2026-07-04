#!/usr/bin/env sage
# merge_results.py
# Merges cm_results_worker*.csv into a single cm_results_merged.csv

import glob, csv, os

pattern = "cm_results_worker*of*.csv"
files = sorted(glob.glob(pattern))

if not files:
    print(f"No files matching '{pattern}' found in {os.getcwd()}")
    print("Looking for files in parent directory...")
    pattern = "../cm_results_worker*of*.csv"
    files = sorted(glob.glob(pattern))

print(f"Found {len(files)} worker files: {files}\n")

all_rows = []
header = None

for fname in files:
    with open(fname) as f:
        reader = csv.reader(f)
        f_header = next(reader)
        if header is None:
            header = f_header
        n = 0
        for row in reader:
            all_rows.append(row)
            n += 1
    print(f"  {fname}: {n} rows")

all_rows.sort(key=lambda r: (r[0], int(r[3])))  # Sort by family, then d

outfile = "cm_results_merged.csv"
with open(outfile, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    w.writerows(all_rows)

print(f"\nMerged {len(all_rows)} rows → {outfile}")

# Summary by family
from collections import Counter
by_family = Counter(r[0] for r in all_rows)
ok_by_family = Counter(r[0] for r in all_rows if r[4] == "OK")
print("\nPer-family summary:")
for fam in sorted(by_family.keys()):
    print(f"  {fam}: {ok_by_family[fam]}/{by_family[fam]} OK")
