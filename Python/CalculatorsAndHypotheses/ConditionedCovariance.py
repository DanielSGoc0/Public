# Calculates conditioned covariance of a fixed coordinates method
import numpy as np
import random

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f} \t'.format}, linewidth=200000, threshold=np.inf)

n = 3
N = 6

X = np.zeros((N + 1, N + 1), dtype='double')
sigma2 = np.zeros((N), dtype='double')

# Indices 0 to n-1 are for some whatever points.
# Indices n to N-1 are for evaluation points.
# Index N is for output

for i in range(n):
	for j in range(n):
		# X[i][j] = 0.0
		X[i][j] = random.random()
for j in range(N):
	X[N][j] = random.random()
for i in range(n, N + 1):
	for j in range(i):
		X[i][j] = X[N][j]
		# X[i][j] = 1.0 / np.sqrt(j + 1)
		# X[i][j] = random.random()

print(X)

for i in range(0, N):
	# sigma2[i] = i
	sigma2[i] = random.random()

Sigma = np.zeros((N + 1, N + 1), dtype='double')
for i in range(N + 1):
	for j in range(N + 1):
		Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
for i in range(0, N):
	Sigma[i][i] += sigma2[i]

ANS = Sigma[N][0] - Sigma[N, n:N] @ np.linalg.inv(Sigma[n:N, n:N]) @ Sigma[n:N, 0]




a = np.zeros((N), dtype='double')
for i in range(n, N):
	a[i] = np.exp(-X[N][i]**2 / 2.0)

gamma2 = np.zeros((N), dtype='double')
kappa = np.zeros((N), dtype='double')
gamma2[n] = sigma2[n]
kappa[n] = 0.0

for i in range(n, N - 1):
	gamma2[i + 1] = sigma2[i + 1] * (1 - a[i]**2 * sigma2[i + 1] / (1.0 + gamma2[i] + a[i]**2 * (sigma2[i + 1] - 1.0)))
	kappa[i + 1] = (1.0 + kappa[i] - a[i]**2) * sigma2[i + 1] / (1.0 + gamma2[i] + a[i]**2 * (sigma2[i + 1] - 1.0))

print(ANS)
print(Sigma[0][N] * (gamma2[N - 1] - kappa[N - 1]) / (1.0 + gamma2[N - 1]))
