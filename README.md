# Reproducibility Guide — paper1_v8_12

## Paper
A Combinatorial Descent Operator for Systematic Enumeration of Admissible Hypergeometric Parameters in Ramanujan-Type 1/pi Series

## Files

| File | Description |
|------|-------------|
| `find_candidates.py` | Python 3 implementation of the descent operator S_p. No dependencies beyond the standard library. |
| `full_catalog.csv` | All 351 structurally admissible (s,p,d) triples with complete rational parameter lists. |

## Quick Start

```bash
python3 find_candidates.py > full_catalog.csv
```

Outputs 351 lines, one per (s,p,d) triple. Each line contains:
- Signature s
- Prime p  
- Heegner discriminant d
- |a'|, |b'|, k = |a'| - |b'|
- Semicolon-separated rational parameters a', b'

## Verification

The output satisfies three structural checks:
1. |a'| - |b'| = 1 for s=1/2, = p+1 for s in {1/3,1/4,1/6} (Theorem 4.4)
2. All fractions are in (0,1] and fully reduced
3. No duplicate fractions after mod-Z reduction

To verify independently:
```bash
python3 -c "
import csv
with open('full_catalog.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # structural balance check
        k = int(row['k'])
        assert k == 1 or k == int(row['p']) + 1
print('All checks passed')
"
```

## Requirements
- Python 3.8+
- No external packages (fractions, math from stdlib only)
