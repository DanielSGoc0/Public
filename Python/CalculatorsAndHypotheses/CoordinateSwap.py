# A simple program
# Calculates Cholesky transform of A,
# as well as Cholesky transform of A',
# where A' has rows/columns i and i+1 swapped.
# We assume A has only 1s on diagonal.
import numpy as np

np.set_printoptions(suppress=True, linewidth=200000, threshold=np.inf)

N = 8
k = 2
A = np.random.random((N, N))
A = A @ A.T

factors = np.zeros((N), dtype='double')
for i in range(N):
	factors[i] = 1.0 / np.sqrt(A[i][i])

for i in range(N):
	for j in range(N):
		A[i][j] *= factors[i] * factors[j]

T = np.eye(N, dtype='double')
T[k][k] = 0.0
T[k+1][k+1] = 0.0
T[k+1][k] = 1.0
T[k][k+1] = 1.0

A2 = T @ A @ T

print(A)
print(A2)

print(np.linalg.cholesky(A))
print(np.linalg.cholesky(A2))

print(np.linalg.cholesky(A) - np.linalg.cholesky(A2))
