# Checks whether various inequalities for Cholesky Decomposition hold

import numpy as np
import random

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)


MIN_v = 10000000
MAX_v = -10000000

N = 4
D = 2
h = 0.001
SCALE = 3.0

worst = None

for k in range(10000):
	X = np.array([[random.random() * SCALE for j in range(D)] for i in range(N + 1)])
	X[N] = X[N - 1]
	X[N][0] += h
	# X[0][0] = 0.0
	# X[0][1] = 0.0
	# X[1][0] = 0.0 

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2)
	try:
		C = np.linalg.cholesky(Sigma)

		w1 = C[N - 1][N - 1]**2
		w2 = C[N][N - 1]**2 + C[N][N]**2

		v = (w1 - w2)/np.sqrt(np.dot(X[N] - X[N - 1], X[N] - X[N - 1]))
		MIN_v = min(MIN_v, v)
		if MAX_v < v:
			MAX_v = max(MAX_v, v)
			worst = X
	except:
		print(k)

print("MAX: ", MAX_v)
print("MIN: ", MIN_v)
print("MAX AT:")
print(worst)
