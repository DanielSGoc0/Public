# ChatGPT generated code.
# Checks if "eigenvectors oscillate" property
# is preserved under multiplication by diagonal
# matrices
import numpy as np
import itertools
import random
from scipy.stats import ortho_group

np.set_printoptions(suppress=True, precision=4)

def sign_changes(vec):
    s = np.sign(vec)
    changes = 0
    for i in range(len(s)-1):
        if s[i]*s[i+1] < 0:
            changes += 1
    return changes

def has_oscillation_property(M):
    w, V = np.linalg.eigh(M)
    if np.min(np.abs(np.diff(w))) < 1e-8:
        return False, w, V
    sc = [sign_changes(V[:,i]) for i in range(V.shape[1])]
    return sc == list(range(M.shape[0]))[::-1], w, V

def loses_property_after_D(M, D):
    M2 = D @ M @ D
    ok, w2, V2 = has_oscillation_property(M2)
    return (not ok), w2, V2

# Search parameters
n = 5
max_attempts = 500000
found = False
np.random.seed(1)

for attempt in range(max_attempts):
    V = ortho_group.rvs(dim=n)
    eigs = sorted([np.random.random() for i in range(n)])
    M = V @ np.diag(eigs) @ V.T
    # Check all entries positive
    if np.all(M > 1e-8):
        # print(M)
        osc, _, _ = has_oscillation_property(M)
        if osc:
            print("HAS OSCILLATION PROPERTY")
            print(M)
            # Try a few random positive diagonals
            for _ in range(500):
                D = np.diag(np.random.uniform(0.2, 2.0, n))
                loses, w2, V2 = loses_property_after_D(M, D)
                if loses:
                    found = True
                    print("Found candidate after", attempt+1, "attempts")
                    print("M (all entries positive):\n", M)
                    print("\nD:\n", D)
                    print("\nEigenvalues of M:", eigs)
                    print("Sign changes of M eigenvectors:", list(range(n)))
                    print("\nEigenvalues of D M D:", np.linalg.eigvalsh(D@M@D))
                    sc_new = [sign_changes(V2[:,i]) for i in range(n)]
                    print("Sign changes of D M D eigenvectors:", sc_new)
                    break
    if found:
        break

if not found:
    print("No counterexample found in", max_attempts, "attempts. Try increasing max_attempts or adjusting parameters.")
