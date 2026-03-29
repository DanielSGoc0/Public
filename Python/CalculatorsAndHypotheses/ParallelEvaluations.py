# Used to check if "normalizing" the evaluation matrix
# always increases the value of <XN, CN>
import numpy as np
import scipy
import random
import scipy.linalg
from pathlib import Path
import itertools
import time

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f} \t'.format}, linewidth=200000, threshold=np.inf)

n = None
l = None
L = None
N = None


# Here is how indices work:
# 0:n corresponds to initial data
# l[i]:l[i + 1] corresponds to parallel evaluations
# Thus l[0] = n and l[L] = N
# N is the output point


def phi(x, y):
	return np.exp(-np.dot(x - y, x - y) / 2.0)
	# return 3.0 * np.exp(-np.dot(x - y, x - y)) + np.exp(-np.dot(x - y, x - y) / 4.0)


def calc(X, sigma2):
	global n, l, L, N

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(X[i], X[j])
		Sigma[i][i] += sigma2[i]
	C = scipy.linalg.cholesky(Sigma, lower=True)

	RES = 0.0
	for i in range(L):
		RES += np.sqrt(np.dot(C[N, l[i]:l[i + 1]], C[N, l[i]:l[i + 1]]) * np.dot(X[N, l[i]:l[i + 1]], X[N, l[i]:l[i + 1]]))
	return RES


def normalize(X, sigma2):
	global n, l, L, N

	X2 = np.zeros((N + 1, N + 1), dtype='double')
	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	C = np.zeros((N + 1, N + 1), dtype='double')

	for i in range(N + 1):
		for j in range(n):
			X2[i][j] = X[i][j]

	for i in range(n):
		for j in range(n):
			Sigma[i][j] = phi(X[i], X[j])
		Sigma[i][i] += sigma2[i]
	C[:n, :n] = scipy.linalg.cholesky(Sigma[:n, :n], lower=True)

	for j in range(n):
		Sigma[N][j] = phi(X[N], X[j])
		C[N][j] = (Sigma[N][j] - np.dot(C[N, :j], C[j, :j])) / C[j][j]

	last_s = 100000000.0
	for k in range(L):
		for i in range(l[k], l[k + 1]):
			for j in range(l[k + 1]):
				Sigma[i][j] = phi(X2[i], X2[j])
				Sigma[j][i] = Sigma[i][j]
			Sigma[i][i] += sigma2[i]
			Sigma[N][i] = phi(X2[N], X2[i])
			Sigma[i][N] = Sigma[N][i]
		
		for i in range(l[k], l[k + 1]):
			for j in range(i):
				C[i][j] = (Sigma[i][j] - np.dot(C[i, :j], C[j, :j])) / C[j][j]
			C[i][i] = np.sqrt(max(sigma2[i], Sigma[i][i] - np.dot(C[i, :i], C[i, :i])))
			C[N][i] = (Sigma[N][i] - np.dot(C[N, :i], C[i, :i])) / C[i][i]

		# Finally, we set X2 to be parallel to C.
		for i in range(l[k + 1], N + 1):
			s = np.sqrt(np.dot(X[i, l[k]:l[k + 1]], X[i, l[k]:l[k + 1]]) / np.dot(C[N, l[k]:l[k + 1]], C[N, l[k]:l[k + 1]]))
			X2[i, l[k]:l[k + 1]] = C[N, l[k]:l[k + 1]] * s

			if i == N:
				if s < last_s:
					last_s = s
				else:
					return(X, False)

	return (X2, True)



# for tries in range(1000):
# 	n = 0
# 	L = 2
# 	l = [0, 2, 4]
# 	N = 4

# 	X = np.zeros((N + 1, N + 1), dtype='double')
# 	for i in range(N + 1):
# 		for j in range(n):
# 			X[i][j] = 1.0 - 2.0 * random.random()
# 			# X[i][j] = random.random()
# 	for k in range(L):
# 		for i in range(l[k + 1], N + 1):
# 			for j in range(l[k], l[k + 1]):
# 				X[i][j] = 1.0 - 2.0 * random.random()
# 				# X[i][j] = random.random()
	
# 	sigma2 = np.zeros((N + 1), dtype='double')
# 	for i in range(N):
# 		sigma2[i] = random.random() / (1.0 - random.random())

# 	ANS1 = calc(X, sigma2)
# 	X2 = normalize(X, sigma2)
# 	ANS2 = calc(X2, sigma2)

# 	EPS = 0.000001

# 	if ANS2 + EPS < ANS1:
# 		print("SHIT")
# 		print(ANS1, ANS2)
# 		print(X)
# 		print(X2)

# 		# Sigma = np.zeros((N + 1, N + 1), dtype='double')
# 		# for i in range(N + 1):
# 		# 	for j in range(N + 1):
# 		# 		Sigma[i][j] = phi(X[i], X[j])
# 		# 	Sigma[i][i] += sigma2[i]
# 		# C = scipy.linalg.cholesky(Sigma, lower=True)



# 		break
# 	else:
# 		print(tries)
# 		print(ANS1, ANS2)



