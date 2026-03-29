import numpy as np
import scipy
import math

def f(w, N, V, X):
	dw = np.zeros(N, dtype='double')
	for i in range(N):
		dw[i] = math.exp(np.dot(X[i], w - X[i]/2.0))
	return scipy.linalg.solve_triangular(V, dw, lower=True)

def error(X, S):
	X = np.array(X)
	N = X.shape[0]

	M = np.zeros((N, N), dtype='double')
	for i in range(N):
		for j in range(N):
			M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	V = np.linalg.cholesky(M)
	# print(V)
	# print(scipy.linalg.solve_triangular(V, np.identity(N), lower=True))

	sol = scipy.integrate.solve_ivp(lambda t, w: f(w, N, V, X), (0, S), np.zeros((N), dtype='double'), method='DOP853')

	RES = -2*math.log(1 - S*S/2)
	# RES = 0.0
	for i in range(N):
		RES -= sol.y[i][-1]**2
		# print(sol.y[i][-1])
	return RES

def get_X(N, v):
	X = np.zeros((N, N), dtype='double')
	k = 0
	for i in range(N):
		for j in range(i):
			X[i][j] = v[k]
			k += 1
	return X

def minimize(N, S):
	v0 = []
	for i in range(N - 1):
		for _ in range(i):
			v0.append(0.0)
		v0.append(1.0)

	fun = lambda v: error(get_X(N, v), S)

	print(scipy.optimize.minimize(fun, v0))

# v = [-1.0, 1.0, 0.01]
# print(error(get_X(3, v), 0.1))
# print(error(get_X(3, v), 0.2))
# print(error(get_X(3, v), 0.3))
# print(error(get_X(3, v), 0.4))

def get_S(T):
	return math.sqrt(2) * math.tanh(T * math.sqrt(2) / 2.0)

print(error([[0.0, 0.0, 0.0], [0.4385, 0.0, 0.0], [0.6541, 0.2786, 0.0]], get_S(1)))
# print(error(get_X(5, [0.2641, 0.5078, 0.1668, 0.7518, 0.3062, 0.1462, 0.9371, 0.4536, 0.2039, 0.1123]), 1))
# print(error(get_X(3, [0.5229, 0.7670, 0.3853]), 1))

# minimize(3, get_S(1))
