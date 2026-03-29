# Calculates the optimal convergence rate

import numpy as np

alpha = 0.9
N = 100

a = np.zeros(N + 1)
a[0] = np.exp(-alpha**2 / 2.0) * alpha

S = 0.0
M = alpha**2

for n in range(N):
	print(a[n])
	
	S += a[n]**2
	a[n + 1] = alpha * np.exp((S - M) / 2.0) * np.sqrt(1.0 - np.exp(-a[n]**2))

print(S)
print(M)
