# Calculates the average determinant of N unit vectors chosen uniformly at random
# from unit sphere in M dimensions. The determinant is considered to be taken from
# the projection of vectors onto their common subspace.

# The problem has been solved completely, and the average det is given by the formula
# val(N, M) below.

import numpy as np
import math
from scipy.special import gamma

# N is the amount of vectors
# M is the dimension of space

def random_unit_vectors(N, M):
	V = np.random.normal(0, 1, (N, M))
	for i in range(N):
		V[i] = V[i] / math.sqrt(np.dot(V[i], V[i]))
	return V

def random_det(N, M):
	V = random_unit_vectors(N, M)
	A = np.inner(V, V)
	return math.sqrt(np.linalg.det(A))

def val(N, M):
	return gamma(M/2.0) / gamma((M - N + 1)/2.0) * (gamma(M/2.0) / gamma((M + 1)/2.0))**(N - 1)

N = 2
M = 3
SUM = 0.0
COUNT = 0.0
PROPOSED = val(N, M)

for i in range(100000000):
	SUM += random_det(N, M)
	COUNT += 1
	print(PROPOSED)
	print(i, '\t', SUM / COUNT)

# 1 vector in R1:
# 1.0

# 1 vector in R2:
# 1.0

# 2 vectors in R2:
# 2/pi

# 1 vector in R3:
# 1.0

# 2 vectors in R3:
# pi/4?

# 3 vectors in R3:
# pi/8?
