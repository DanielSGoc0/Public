# Covariance is exp(-||x - y||^2)
# Implementation of searching a fixed point.

import numpy as np
import scipy
import math
import random

def get_v(N, X):
	v = []
	for i in range(N + 1):
		for j in range(i):
			v.append(X[i][j])
	return v


def get_X(N, v):
	X = np.zeros((N + 1, N + 1), dtype='double')
	k = 0
	for i in range(N + 1):
		for j in range(i):
			X[i][j] = v[k]
			k += 1

	# X[1][0] = 1.0
	# X[2][0] = 1.3 
	# X[2][1] = 1.7 
	# X[3][0] = 0.5 
	# X[3][1] = 0.8 
	# X[3][2] = 1.7 
	return X


# ===============================================================================

def get_sum(X, I, s, t):
	SUM = 0.0
	for j in range(s):
		for k in range(t):
			SUM += X[s][j] * X[t][k] * np.exp(-I[j][j] - I[k][k] + 2*I[j][k])
	return SUM

# X is strictly lower diagonal
def calc(X, alpha):
	X = np.array(X)
	N = X.shape[0] - 1

	I = np.zeros((N + 1, N + 1))

	for s in range(1, N + 1):
		for t in range(1, s + 1):
			SUM = 0.0

			for j in range(s):
				for k in range(t):
					SUM += X[s][j] * X[t][k] * np.exp(-I[j][j] - I[k][k] + 2*I[j][k])

			I[s][t] = SUM
			I[t][s] = SUM

	X[N][N] = -alpha
	ANS = 0.0
	for j in range(N + 1):
		for k in range(N + 1):
			ANS += X[N][j] * X[N][k] * np.exp(-I[j][j] - I[k][k] + 2*I[j][k])

	print(ANS)
	print(X)
	print(I)
	# print(X[1][0] + (2*X[1][0]*X[2][0] - 1)*alpha*np.exp(-I[2][2]))
	# print(-X[2][0] - 2*alpha*X[2][0]**2 * np.exp(-I[2][2]) - 2*alpha*X[2][1]*(X[2][0] - X[1][0]) * np.exp(-I[1][1] - I[2][2] + 2*I[1][2]))
	# print((X[2][0] - X[1][0]) + X[2][1] * np.exp(-X[1][0]*X[1][0]))
	# print(X[2][0] * np.exp(-I[2][2]) + 2*X[1][0]*X[2][1]*(X[2][0] - X[1][0])*np.exp(-I[1][1]-I[2][2]+2*I[1][2]))
	# print(X[2][1] * np.exp(-I[2][2]) + X[1][0] * np.exp(-I[1][1]-I[2][2]+2*I[1][2]) * (-1 + 2*X[2][1]**2))
	return ANS


# =========================================================================================

def initial_guess(N, alpha):
	X = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N):
		# X[i + 1][i] = alpha
		X[i + 1][i] = random.random()

	return X


def minimize(N, alpha):
	v0 = initial_guess(N, alpha)
	print(v0)

	fun = lambda v: calc(get_X(N, v), alpha)

	RES = scipy.optimize.minimize(fun, get_v(N, v0), options={"gtol": 0.00000001})
	print("=================================================================")
	print(RES)
	print(get_X(N, RES.x))
	fun(RES.x)


# print(calc([[0.0, 0.0, 0.0], [0.4385, 0.0, 0.0], [0.6541, 0.2786, 0.0]], 0.5))

minimize(5, 1.0)

# v0 = [0.5895, 0.8148, 0.6386]
# calc(get_X(3, v0), 30, 0.01)
