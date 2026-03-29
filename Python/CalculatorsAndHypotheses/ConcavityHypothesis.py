# Used to check if <C_N, X_N> is concave in X.
import numpy as np
import scipy
import random
import scipy.linalg
from pathlib import Path
import itertools
import time

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f} \t'.format}, linewidth=200000, threshold=np.inf)

n = None
m = None
N = None


# Here is how indices work:
# 0:n corresponds to initial data
# n:m corresponds to first parallel evaluation
# m:N corresponds to latter evaluations
# N is the output point

def phi(x, y):
	r2 = np.dot(x - y, x - y)
	return np.exp(-r2 / 2.0)

def get_Sigma(X, sigma2):
	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			# Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
			Sigma[i][j] = phi(X[i], X[j])

	for i in range(n):
		Sigma[i][i] += sigma2[i]
	for i in range(m, N):
		Sigma[i][i] += sigma2[i]

	return Sigma


# BUT X on input has shape [0:(N + 1), 0:(N + 1 - (m - n - 1))],
# as X[0:(N + 1), n] is expanded into X[0:(n + 1), n:m]
def calc(X, Sigma, sigma2):
	Sigma_copy = np.copy(Sigma)
	for i in range(n, m):
		Sigma_copy[i][i] += sigma2[i]
	C = scipy.linalg.cholesky(Sigma_copy, lower=True)

	return np.sqrt(np.dot(C[N, n:m], C[N, n:m])) * np.abs(X[N][n]) + np.dot(np.abs(C[N, m:N]), np.abs(X[N, (n+1):(N - (m - n - 1))]))
	# return np.sqrt(np.dot(C[N, n:m], C[N, n:m])) * X[N][n]



for tries in range(1000):
	n = 3
	m = 5
	N = 6

	X = np.zeros((N + 1, N + 1 - (m - n - 1)), dtype='double')
	sigma2 = np.zeros((N), dtype='double')

	for i in range(N + 1):
		for j in range(n):
			# X[i][j] = 1.0 - 2.0 * random.random()
			X[i][j] = random.random()
	for i in range(m, N + 1):
		for j in range(i - (m - n - 1)):
			X[i][j] = random.random()
	
	P = 1.0

	for i in range(N):
		sigma2[i] = random.random() / (1.0 - random.random())

	sigma2_A = np.copy(sigma2)
	sigma2_B = np.copy(sigma2)
	sigma2_C = np.copy(sigma2)
	for i in range(n, m):
		sigma2_C[i] = random.random() / (1.0 - random.random())

	sum_A = 0.0
	sum_C = 0.0
	for i in range(n, m):
		sum_A += 1.0 / sigma2_A[i]
		sum_C += 1.0 / sigma2_C[i]
	for i in range(n, m):
		sigma2_A[i] *= sum_A / P
		sigma2_C[i] *= sum_C / P
		sigma2_B[i] = 2.0 / (1.0 / sigma2_A[i] + 1.0 / sigma2_C[i])

	# print(sigma2_A)
	# print(sigma2_B)
	# print(sigma2_C)

	Sigma = get_Sigma(X, sigma2)

	ANS_A = calc(X, Sigma, sigma2_A)
	ANS_B = calc(X, Sigma, sigma2_B)
	ANS_C = calc(X, Sigma, sigma2_C)

	if ANS_B < min(ANS_A, ANS_C):
	# if ANS_A + ANS_C < 2.0 * ANS_B:
		print("SHIT")
		print(ANS_A, ANS_B, ANS_C)
		break
	else:
		print(tries)
		print(ANS_A, ANS_B, ANS_C)



