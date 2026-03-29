# Covariance is exp(-||x - y||^2 / 2)
# A certain hypothesis test

import numpy as np
import scipy
import random

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

err2 = 10.0

def get_v(N, X, modifiable):
	v = []
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				v.append(X[i][j])
	return v


def get_X(N, v, X0, modifiable):
	X = np.zeros((N + 1, N + 1), dtype='double')
	k = 0
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				X[i][j] = v[k]
				k += 1
			else:
				X[i][j] = X0[i][j]

	return X


# ===============================================================================

# X is strictly lower diagonal
def calc(X):
	global err2

	N = len(X) - 1
	X = np.array(X)

	M = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			M[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		M[i][i] += err2
	Chol = scipy.linalg.cholesky(M, lower=True)

	return (np.dot(X[N], Chol[N]), np.dot(X[N - 1], Chol[N - 1]), Chol[N - 1][N - 1]**2)


# =========================================================================================

def initial_guess(N):
	X = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i):
			X[i][j] = random.random()

	return X


def minimize(x0):
	N = 3
	h = 0.00001

	X0 = initial_guess(N)
	X0[1][0] = x0

	modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
	for i in range(2, N + 1):
		for j in range(i):
			modifiable[i][j] = True

	fun = lambda v: -calc(get_X(N, v, X0, modifiable))[0]

	RES = scipy.optimize.minimize(fun, get_v(N, X0, modifiable), options={"gtol": 0.00000000000001})
	X0 = get_X(N, RES.x, X0, modifiable)

	C = [0.0, 0.0, 0.0]
	F = [0.0, 0.0, 0.0]

	vals = calc(X0)
	F[0] = vals[1]
	C[0] = vals[2]

	X0[2][0] += h
	vals = calc(X0)
	F[1] = vals[1]
	C[1] = vals[2]
	X0[2][0] -= h

	X0[2][1] += h
	vals = calc(X0)
	F[2] = vals[1]
	C[2] = vals[2]
	X0[2][1] -= h

	print(x0, (F[1] - F[0]) / (C[1] - C[0]), (F[2] - F[0]) / (C[2] - C[0]))
	# print((F[1] - F[0]) / h)
	# print((C[1] - C[0]) / h)
	# print((F[2] - F[0]) / h)
	# print((C[2] - C[0]) / h)

for k in range(1, 200):
	minimize(k / 100.0)
