# Calculates a random determinant of matrix as given below.
# Fairly straightforward

import numpy as np
import random
import scipy
import itertools

import scipy.special

N = 3
P = 3
D = 2

def random_determinant():
	global N, P, D

	X = np.array([[random.normalvariate(0.0, 1.0) for d in range(D)] for n in range(N)])

	A = np.zeros((N, N), dtype='double')

	for i in range(N):
		for j in range(N):
			A[i][j] = np.dot(X[i], X[j])**P

	return np.linalg.det(A)

for k in range(100):
	print(random_determinant())
	