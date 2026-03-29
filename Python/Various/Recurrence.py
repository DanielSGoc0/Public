# Calculates optimal coefficients for TrueGradient.
# Turns out they are not numerically stable

import numpy as np

N = 100

a = [0.0 for i in range(N)]
c = [0.0 for i in range(N)]
w = [0.0 for i in range(N)]

a[0] = 1.28650116
# a[0] = 1.3
c[0] = a[0]
w[0] = np.exp(-a[0]**2/2.0)
print(0, a[0], c[0], w[0])


c[1] = (c[0] - 1/c[0])/w[0]
a[1] = c[1] * np.sqrt(1 - w[0]**2)
w[1] = np.exp(-c[1]**2 * (1 - w[0]**2) / 2.0)
print(1, a[1], c[1], w[1])


for k in range(2, N):
	c[k] = (w[k - 2]/c[k - 2] + c[k - 1] - 1.0/c[k - 1])/w[k - 1]
	a[k] = c[k] * np.sqrt(1 - w[k - 1]**2)
	w[k] = np.exp(-c[k]**2 * (1 - w[k - 1]**2) / 2.0)
	print(k, a[k], c[k], w[k])
