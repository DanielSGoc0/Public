# Checks whether various inequalities for Cholesky Decomposition hold

import numpy as np

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)


n = 2
d = 1
sigma = 2.5
c = 1.0

for k in range(100000):
	if k % 1000 == 0:
		print(k)
	X = np.array([[np.random.random() * sigma for j in range(d)] for i in range(n + 2)])

	Sigma = np.zeros((n + 2, n + 2), dtype='double')
	for i in range(n + 2):
		for j in range(n + 2):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2)
	for i in range(n):
		p = np.random.random()
		Sigma[i][i] += p/(1.0 - p) * c
	C = np.linalg.cholesky(Sigma)

	cov = C[n+1][n] * C[n][n]
	if cov > Sigma[n+1][n]:
		print(X)
		print(cov)
		print(Sigma[n+1][n])
		break
	elif cov < -Sigma[n+1][n]:
		print(k, '\t.')



