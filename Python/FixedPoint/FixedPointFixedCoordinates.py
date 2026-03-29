# Covariance is exp(-||x - y||^2 / 2)
# Implementation of searching for a fixed point problem.
# Adjusted for one of the optimality results! (i.e. fixed coordinates method)

import numpy as np
import scipy
import random

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=130)

# ===============================================================================
err2 = 1.0

def calc(X, alpha):
	global err2

	N = len(X)

	mags = np.zeros((N + 1))
	for i in range(N):
		mags[i + 1] = mags[i] + X[i]**2

	M = np.zeros((N + 1, N + 1))
	for i in range(N + 1):
		for j in range(N + 1):
			M[i][j] = np.exp((mags[min(i, j)] - mags[max(i, j)]) / 2.0)
		M[i][i] += err2
	Chol = scipy.linalg.cholesky(M, lower=True)

	ANS = alpha**2 * Chol[N][N]**2
	# ANS = 0
	for i in range(N):
		ANS += (X[i] - alpha * Chol[N][i])**2

	# print(Chol)
	print(ANS - err2 * alpha**2)
	print(X)
	return ANS


# =========================================================================================

def initial_guess(N, alpha):
	X = np.zeros((N), dtype='double')
	for i in range(N):
		X[i] = random.random()

	return X


def minimize(N, alpha):
	X0 = initial_guess(N, alpha)
	print(X0)

	fun = lambda X: calc(X, alpha)

	RES = scipy.optimize.minimize(fun, X0, options={"gtol": 0.00000001})
	print("=================================================================")
	print(RES)
	print(RES.x)
	fun(RES.x)

minimize(5, 1.0)

# calc([0.5, 0.6, 0.7, 0.9, 0.3], 1.0)
# calc([0.5], 1.0)
