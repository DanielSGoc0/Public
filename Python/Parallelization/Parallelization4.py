# Covariance is exp(-||x - y||^2 / 2)
# This one is for checking the convexity of conditioned covariance.
# PROVEN!

import numpy as np
import scipy
import random

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)



def calc(X, precision):
	X = np.array(X)
	N = len(X) - 1

	# The first N elements are evaluation points.
	# The N-th point is the target point.

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		Sigma[i][i] += 1.0/precision[i]
		
	return Sigma[N][N] - np.dot(Sigma[N][:N], scipy.linalg.solve(Sigma[:N, :N], Sigma[N][:N], assume_a='pos'))

def generate(N, D, TOTAL):

	SUM = 0.0
	precision = np.zeros((N + 1), dtype='double')
	for i in range(N):
		precision[i] = np.exp(3.0 - 6.0*random.random())
		SUM += precision[i]
	for i in range(N):
		precision[i] *= TOTAL/SUM
	precision[N] = 1.0

	X = np.array([[2.0 - random.random()*4.0 for j in range(D)] for i in range(N + 1)])
	for j in range(D):
		X[N][j] = 0.0

	return (X, precision)


N = 2
D = 1

for i in range(1000):
	(X0, precision0) = generate(N, D, 1.0)
	(_, precision1) = generate(N, D, 1.0)

	precision2 = (precision0 + precision1) / 2.0

	val0 = calc(X0, precision0)
	val1 = calc(X0, precision1)
	val2 = calc(X0, precision2)

	if (val0 + val1)/2.0 < val2:
		print(X0)
		print(precision0)
		print(precision1)
		print(precision2)
		break

	print((val0 + val1)/2.0 - val2)
