# This is a very important program.
# We analyze the behavior of weights w_k under
# the optimality conditions applied to fixed
# coordinates method (this is the errored case)
import numpy as np
from scipy.linalg import solve_triangular

# Indices 0:N are for evaluations. N is for output
N = 1000
R2 = 0.5
rho2 = 0.5
p = np.array([0.01 for i in range(N)])

# h2 is defined in such a way that X[N][k] = h * C[N][k]
h2 = R2 / (1.0 - np.exp(-R2) * rho2)
print("h2 must be in interval [0, 1]:", h2)



# w2[k] = C[k][k]**2 - sigma2[k]
# x2[k] = X[N][k]**2
# r2[k] = x2[k] + x2[k + 1] + ... + x2[N - 1],
# where informally x2[N - 1] is such that r2[0] = R2
x2 = np.zeros((N - 1), dtype='double')
w2 = np.zeros((N), dtype='double')
r2 = np.zeros((N), dtype='double')

w2[0] = 1.0 - rho2
r2[0] = R2

for k in range(N - 1):
	x2[k] = h2 * np.exp(-r2[k]) * w2[k]**2 / (w2[k] + 1.0 / p[k])
	w2[k + 1] = 1.0 - np.exp(-x2[k]) * (1.0 - w2[k] / (1.0 + p[k] * w2[k]))
	r2[k + 1] = r2[k] - x2[k]


# Next calculate coefficients gamma2
gamma2 = np.zeros((N), dtype='double')
gamma2[0] = 1.0 / p[0]

for k in range(N - 1):
	a2 = np.exp(-x2[k])
	gamma2[k + 1] = (1.0 + gamma2[k] - a2) / (a2 + p[k + 1] * (1.0 + gamma2[k] - a2))


# and finally, calculate coefficients
kappa = np.zeros((N, N), dtype='double')

for i in range(N):
	kappa[i][i] = gamma2[i]

	for j in range(i, N - 1):
		a2 = np.exp(-x2[j])
		kappa[i][j + 1] = kappa[i][j] / (a2 + p[k + 1] * (1.0 + gamma2[k] - a2))



# Our assumption for weights w[k] is that:
# m = sum_{j = i}^{N - 1} w[j] * (kappa[i][j] / (1 + gamma2[j]))**2
M = np.zeros((N, N), dtype='double')
for i in range(N):
	for j in range(i, N):
		M[i][j] = (kappa[i][j] / (1.0 + gamma2[j]))**2

v = np.zeros((N), dtype='double')
for i in range(N):
	v[i] = np.exp(r2[i])

print(gamma2)

# print(v)

# # With that, we assume m = 1 and obtain w.
# w = solve_triangular(M, v, lower=False)

# print(w)
