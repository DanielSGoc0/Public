# Calculates the posterior variance, where initial data is given
# at a unit circle.

import numpy as np
import scipy
import math
import random

np.set_printoptions(suppress=True)


def Cov(x):
	return math.exp(-np.dot(x, x) / 2.0)

def get_matrix(X, Z):
	X = np.array(X)
	Z = np.array(Z)
	D = X.shape[0]
	K = Z.shape[0]

	Sigma11 = [[Cov(X[i] - X[j]) for j in range(D)] for i in range(D)]
	Sigma12 = [[Cov(X[i] - Z[j]) for j in range(K)] for i in range(D)]
	Sigma21 = [[Cov(Z[i] - X[j]) for j in range(D)] for i in range(K)]
	Sigma22 = [[Cov(Z[i] - Z[j]) for j in range(K)] for i in range(K)]

	SIGMA = Sigma11 - np.dot(Sigma12, scipy.linalg.solve(Sigma22, Sigma21, assume_a = "pos"))
	# SIGMA = np.dot(Sigma12, scipy.linalg.solve(Sigma22, Sigma21, assume_a = "pos"))

	return SIGMA

def equi_dist(K):
	Z = []
	for i in range(K):
		Z.append([math.cos(2*math.pi*i/K), math.sin(2*math.pi*i/K)])
	return Z

def random_Z(K, N):
    Z = np.random.randn(N, K)
    Z /= np.linalg.norm(Z, axis=0)
    return np.transpose(Z)

for k in range(1, 100):
	# SIGMA = get_matrix([[0.0, 0.0], [0.5, 0.0], [1.5, 0.0]], equi_dist(k))
	SIGMA = get_matrix([[0.0, 0.0]], equi_dist(k))
	# SIGMA = get_matrix([[0.0, 0.0], [0.5, 0.0]], random_Z(k))
	# '%.60f'
	print(k)
	# print('%.60f' % SIGMA[0][0])
	print(SIGMA)

# print(random_Z(2, 5))
