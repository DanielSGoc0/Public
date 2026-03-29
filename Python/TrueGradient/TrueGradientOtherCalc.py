# The covariance function is e^(-||x - y||^2 / 2)
# The computation model is different.
# allows for applying error after evaluation.
# Also contains an option to change the formula!
import scipy
import numpy as np
import random

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

err2 = 1000000.0

# =========================================================================================

def get_v(N, X):
	v = []
	for i in range(N + 1):
		for j in range(i):
			v.append(X[i][j])
	return v

def get_v2(N, X):
	v = []
	for j in range(N):
		v.append(X[N][j])
	return v

def get_X(N, v):
	X = np.zeros((N + 1, N + 1), dtype='double')
	k = 0
	for i in range(N + 1):
		for j in range(i):
			X[i][j] = v[k]
			k += 1
	return X

def get_X2(N, v):
	X = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i):
			X[i][j] = v[j]
	return X


# =========================================================================================

# X is strictly lower triangular
def calc(X):
	global err2
	X = np.array(X, dtype='double')
	N = X.shape[0] - 1

	# Sigma_ij = e^(-||Xi - Xj||^2 / 2) + delta_{i, j} err2
	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		Sigma[i][i] += err2

	Chol = scipy.linalg.cholesky(Sigma, lower=True)

	# print(Chol)
	# print(np.linalg.inv(Chol))
	# We calculate f(X_N) only
	RES = np.dot(X[N], Chol[N])
	# RES = np.dot(X[N], Sigma[N])
	# print(np.dot(X[N], Chol[N]) * np.sqrt(err2), np.dot(X[N], Sigma[N]))
	# RES = 0.0
	# for k in range(N + 1):
	# 	RES += np.dot(X[k], Chol[k]) * Chol[k][k]**2
	# RES += np.dot(X[N], Chol[N])
	print(RES)
	return RES

# basically, it's like this: We have a distance squared matrix M.
# Then we consider two matrices: A = Cholesky(M) and B = Cholesky(exp(M))
# The error is then <A[N], B[N]>

def calc_old(X):
	global err2

	X = np.array(X, dtype='double')
	# print("X = ", X)
	N = X.shape[0] - 1

	# Sigma_ij = e^(-||Xi - Xj||^2 / 2)
	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	S = np.zeros((N + 1), dtype='double')

	# Sigma
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		Sigma[i][i] += err2
	Chol = scipy.linalg.cholesky(Sigma, lower=True)

	# S
	S[0] = np.sqrt(1 + err2)
	print(S[0], Chol[0][0])
	for n in range(1, N + 1):
		v = np.array(Sigma[n, :n])
		Sigma_partial = np.array(Sigma[:n, :n])
		S[n] = np.sqrt(1.0 + err2 - np.dot(v, scipy.linalg.solve(Sigma_partial, v, assume_a='pos')))
		print(S[n], Chol[n][n])

	# We calculate f(X_N) only
	F = 0.0

	for n in range(N):
		if n > 0:
			u = np.array(Sigma[N, :n])
			v = np.array(Sigma[n, :n])
			Sigma_partial = np.array(Sigma[:n, :n])
			w = np.dot(u, scipy.linalg.solve(Sigma_partial, v, assume_a='pos'))
			F -= w * X[N][n] / S[n]
			# F -= np.dot(u, scipy.linalg.solve(Sigma_partial, v, assume_a='pos')) * X[N][n]

		F += Sigma[N][n] * X[N][n] / S[n]

	print("F = ", F)
	return F

# A has elements indexed from 0 to N - 1.
def calc_const(A):
	A = np.array(A, dtype='double')
	N = A.shape[0]

	w = np.zeros((N + 1), dtype='double')
	F = np.zeros((N + 1), dtype='double')

	# initial values
	w[0] = 1.0
	F[0] = 0.0

	for n in range(N):
		# calculate w[n + 1] and F[n + 1]
		w[n + 1] = 1 - np.exp(-A[n]**2) * (1 - err2 * w[n] / (err2 + w[n]))
		F[n + 1] = np.exp(-A[n]**2 / 2.0) * (F[n] + A[n] * w[n] / np.sqrt(err2 + w[n]))

	return F[N]


# =========================================================================================

def initial_guess(N):
	X0 = np.zeros((N + 1, N + 1), dtype='double')

	for i in range(N + 1):
		for j in range(N + 1):
			if j < i:
				X0[i][j] = random.random()
			else: 
				X0[i][j] = 0.0

	return X0


def minimize(N):
	X0 = initial_guess(N)
	print(X0)

	fun = lambda v: -calc(get_X(N, v))

	# RES = scipy.optimize.minimize(fun, get_v(N, X0))
	RES = scipy.optimize.minimize(fun, get_v(N, X0), options={"gtol": 0.000000000000000000001})
	
	print("==========================================================")
	print(RES)

	X = get_X(N, RES.x)
	print(X)


def minimize2(N):
	v0 = np.ones((N), dtype='double')

	fun = lambda v: -calc(get_X2(N, v))

	# RES = scipy.optimize.minimize(fun, get_v(N, X0))
	RES = scipy.optimize.minimize(fun, v0, options={"gtol": 0.000000000000000000001})
	
	print("==========================================================")
	print(RES)

	X = get_X2(N, RES.x)
	print(X[-1])


minimize(5)
# minimize2()
# calc([[0, 0, 0, 0], [0.5, 0, 0, 0], [0.3, 0.8, 0, 0], [0.7, -0.2, 0.4, 0]])
# calc_old([[0.00000000, 0.00000000, 0.00000000, 0.00000000], [1.00000000, 0.00000000, 0.00000000, 0.00000000], [0.60000000, 1.58149906, 0.00000000, 0.00000000], [1.47141105, 0.28071358, 0.79723424, 0.00000000]])


# calc_old([[0, 0, 0, 0, 0], [1.1, 0, 0, 0, 0], [1.1, 0.9, 0, 0, 0], [1.1, 0.9, 0.6, 0, 0], [1.1, 0.9, 0.6, 0.4, 0]])
# [[0 0 0 0 0] [1.1 0 0 0 0] [1.1 0.9 0 0 0] [1.1 0.9 0.6 0 0] [1.1 0.9 0.6 0.4 0]]
# print(calc_const([1.1, 0.9, 0.6, 0.4]))
