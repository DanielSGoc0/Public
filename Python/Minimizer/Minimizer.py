# A practical (and naive) implementation of Minimizer resulting from TrueGradient

import numpy as np
import random
from scipy.stats import ortho_group
from scipy.linalg import solve_triangular
import math

err2 = 0
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.4528, 0.0, 0.0, 0.0, 0.0], [0.6651, 0.4350, 0.0, 0.0, 0.0], [0.8443, 0.6241, 0.5030, 0.0, 0.0], [0.9524, 0.7839, 0.6929, 0.6223, 0.0], [1.0, 0.88066503, 0.85910537, 0.86746275, 0.84539765]])
K = X.shape[1]
C = np.zeros((K + 1, K + 1))
h = 1

N = 1
A = np.zeros((N, N), dtype='double')

def f(x):
	global A
	# v = np.dot(A, x)
	# return np.dot(v, v)
	return np.array([math.sin(x[0])])

def df(x):
	global A
	# return np.dot(np.transpose(A), np.dot(A, x))
	return np.array([math.cos(x[0])])

def create_random_matrix():
	global A, N
	A = ortho_group.rvs(dim = N)
	for i in range(N):
		c = random.random()
		c = c/(1 - c)
		for j in range(N):
			A[i][j] *= c

def get_matrix_C():
	global X, C, K

	M = np.zeros((K, K), dtype='double')
	for i in range(K):
		for j in range(K):
			M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]))
		M[i][i] += err2
	V = np.linalg.cholesky(M)
	# print(V)

	C = solve_triangular(np.transpose(V), np.transpose(X), lower=False)
	C = np.transpose(C)
	# print(C)

def iteration(z0):
	global K, C, N, h
	z = np.zeros((K + 1, N))
	z[0] = z0
	ngrads = np.zeros((N, K))
	print(f(z0))

	for k in range(K):
		grad = df(z[k])
		for i in range(N):
			ngrads[i][k] = -grad[i]
		z[k + 1] = z0 + h * np.dot(ngrads, C[k + 1])
		print(f(z[k + 1]))

	return z[-1]


# create_random_matrix()
# print(A)
# A = np.identity(N)
# for i in range(N):
# 	A[i][i] = i
get_matrix_C()
z0 = np.zeros(N)
z0[0] = 1
for i in range(10):
	z0 = iteration(z0)


# x0 = np.zeros(N)
# x0[0] = 1
# x1 = np.zeros(N)
# x1[1] = 1
# print(f(x0))
# print(f(x1))
# print(f(x0 + x1))
