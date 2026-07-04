"""Test modular polynomial computation in SageMath"""
from sage.all import *

# Method 1: Check if PARI/GP has polmodular
try:
    phi_pari = pari.polmodular(11)
    print("PARI polmodular(11) succeeded")
    # phi_pari is a vector of polynomials, need to reconstruct the bivariate polynomial
    print(type(phi_pari))
except Exception as e:
    print(f"PARI method failed: {e}")

# Method 2: Check the modular polynomial database
from sage.databases.db_modular_polynomials import ClassicalModularPolynomialDatabase
import os

# Check what data files exist
sage_root = os.environ.get('SAGE_ROOT', '/mnt/d/miniforge3')
for root, dirs, files in os.walk(sage_root):
    for f in files:
        if 'modular' in f.lower() or 'kohel' in f.lower():
            fp = os.path.join(root, f)
            try:
                sz = os.path.getsize(fp)
                print(f"  {fp} ({sz} bytes)")
            except:
                pass
    # limit depth
    if root.count(os.sep) - sage_root.count(os.sep) > 4:
        dirs.clear()

# Method 3: Check if we can install the database
print("\nChecking for sage-package command...")
import subprocess
try:
    result = subprocess.run(['sage-package', 'list', '--installed'], 
                          capture_output=True, text=True, timeout=30)
    print(result.stdout[:500])
except:
    print("sage-package not available")
