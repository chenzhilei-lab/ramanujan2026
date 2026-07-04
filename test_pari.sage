"""Quick test of PARI polmodular"""
from sage.all import pari
phi = pari.polmodular(11)
print("PARI phi_11 computed successfully")
print(f"Type: {type(phi)}")
print(f"Length: {len(phi)}")
# phi is returned as a vector of polynomials
# For classical modular polynomial, phi = [F, G] where F(X) = scalar * Phi_N(X, Y=Y)
# Actually in latest PARI, polmodular(N) returns a vector [F, G] where
# F(X) = Phi_N(X, Y) and G contains information about the genus
print(f"phi[0] = {phi[0]}")
