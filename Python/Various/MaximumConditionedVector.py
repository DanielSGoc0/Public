# Checks what the largest possible modulus of
# Sigma^{-1} v can be.

import numpy as np
import scipy
import math
import random

np.set_printoptions(suppress=True)


def Cov(x):
	return math.exp(-np.sqrt(np.dot(x, x)) / 2.0)
	# return math.exp(-np.dot(x, x) / 2.0)

def calc():
	N = 1
	D = 1

	X = np.random.random((N + 2, D))
	X[N] = np.zeros((D), dtype='double')
	X[N + 1] = np.ones((D), dtype='double')

	Sigma = np.array([[Cov(X[i] - X[j]) for j in range(N)] for i in range(N)])
	v1 = np.array([Cov(X[N] - X[j]) for j in range(N)])
	v2 = np.array([Cov(X[N + 1] - X[j]) for j in range(N)])

	ANS = Cov(X[N] - X[N + 1]) - v1 @ scipy.linalg.solve(Sigma, v2, assume_a = "pos")
	return (X, ANS)


MAX = 0.0
best_X = None
for k in range(10000):
	(X, val) = calc()
	if abs(val) > MAX:
		MAX = abs(val)
		best_X = X
	print(val, '\t', MAX)

print(best_X)
