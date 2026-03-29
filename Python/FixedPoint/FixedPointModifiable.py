# Covariance is exp(-||x - y||^2 / 2)
# Implementation of searching for a fixed point problem.
# Allows for optimization over a chosen set of method coefficients

import numpy as np
import scipy
import random

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

err2 = 1.0

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
def calc(X, alpha):
	global err2

	N = len(X) - 1
	X = np.array(X)
	# print(X)

	M = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			M[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		M[i][i] += err2
	Chol = scipy.linalg.cholesky(M, lower=True)
	
	ANS = alpha**2 * Chol[N][N]**2
	for i in range(N):
		ANS += (X[N][i] - alpha * Chol[N][i])**2

	# print(Chol)
	# print(X)
	print(ANS)
	return ANS


# =========================================================================================

def initial_guess(N, alpha):
	X = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i):
			X[i][j] = random.random()

	return X


def minimize(N, alpha):
	X0 = initial_guess(N, alpha)
	print(X0)

	modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
	# for j in range(N - 2):
	# 	modifiable[N - 2][j] = True
	# for j in range(N - 1):
	# 	modifiable[N - 1][j] = True
	for i in range(N + 1):
		for j in range(i):
			modifiable[i][j] = True

	fun = lambda v: calc(get_X(N, v, X0, modifiable), alpha)

	RES = scipy.optimize.minimize(fun, get_v(N, X0, modifiable), options={"gtol": 0.00000001})
	print("=================================================================")
	print(RES)
	print(get_X(N, RES.x, X0, modifiable))


# print(calc([[0.0, 0.0, 0.0], [0.4385, 0.0, 0.0], [0.6541, 0.2786, 0.0]], 0.5))

minimize(5, 1.0)

# v0 = [0.5895, 0.8148, 0.6386]
# calc(get_X(3, v0), 30, 0.01)
