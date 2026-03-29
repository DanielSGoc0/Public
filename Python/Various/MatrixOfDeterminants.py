# Creates a matrix whose entries are determinants. Pretty useless.
import numpy as np
import itertools
import random
from scipy.stats import ortho_group

n = 2

V = ortho_group.rvs(dim=n)
eigs = sorted([10*np.random.random() for i in range(n)])
M = V @ np.diag(eigs) @ V.T

A = np.zeros((n, n), dtype='double')

for i in range(n):
	for j in range(n):
		A[i][j] = np.linalg.det(M[min(i, j):(max(i, j) + 1), min(i, j):(max(i, j) + 1)])

print(M)
print(A)


(evals, evecs) = np.linalg.eigh(A)
print(evals)
print(evecs)
