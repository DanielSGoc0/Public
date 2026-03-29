# Needed to check some bug. Turns out everything is alright!

import numpy as np
import random 

import scipy.linalg

N = 3
D = 1

def Cholesky(X):
	global N, error2

	Sigma = np.zeros((N, N), dtype='double')
	for i in range(N):
		for j in range(N):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		Sigma[i][i] += error2[i]
		# Sigma[i][i] += 1.0

	# C = np.linalg.cholesky(Sigma)
	# print(X)
	# print(C)
	# print("===============")
	# print(C[N - 1][N - 3]**2 + C[N - 1][N - 2]**2)
	v12 = Sigma[N-1][N-1] - np.dot(Sigma[N-1][:(N-1)], scipy.linalg.solve(Sigma[:(N-1), :(N-1)], Sigma[N-1][:(N-1)], assume_a='pos'))
	v1 = Sigma[N-1][N-1] - np.dot(Sigma[N-1][:(N-3)], scipy.linalg.solve(Sigma[:(N-3), :(N-3)], Sigma[N-1][:(N-3)], assume_a='pos'))
	# print(v1 - v12)
	# print(v1, v12)

	# return np.linalg.cholesky(Sigma)
	return v12

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

	SUM = 0.0	
	for i in range(N):
		SUM += precision[i]
	print("SUM:", SUM)

	return Sigma[N][N] - np.dot(Sigma[N][:N], scipy.linalg.solve(Sigma[:N, :N], Sigma[N][:N], assume_a='pos'))




for s in range(1):
# for s in range(1):
	# if s % 1000 == 0:
	print("s =", s)
	error2 = [random.random() for k in range(N)]
	error2 = [0.3736601025327182, 0.6721627208456618, 0.5411650405254994]
	# error2 = [0.1 for k in range(N)]
	# X0 = np.array([[random.random()*2.0 for j in range(D)] for i in range(N)])
	X0 = [[1.49161777], [0.37608213], [1.08462918]]
	X0 = np.array(X0)
	# X0[0][0] = 0.0
	# X0[0][1] = 0.0
	# X0[1][1] = 0.0

	opt0 = Cholesky(X0)

	precision0 = [1.0/s for s in error2]

	print(opt0, calc(X0, precision0))



	if True:
		X1 = np.zeros((N, D), dtype='double')
		for i in range(N):
			for j in range(D):
				X1[i][j] = X0[i][j]
		for j in range(D):
			X1[N - 2][j] = X0[N - 3][j]

		opt1 = Cholesky(X1)
		precision1 = [1.0/s for s in error2]
		precision1[N-2] = 0.000001
		precision1[N-3] = precision0[N-2] + precision0[N-3]
		print(opt1, calc(X0, precision1))

		X2 = np.zeros((N, D), dtype='double')
		for i in range(N):
			for j in range(D):
				X2[i][j] = X0[i][j]
		for j in range(D):
			X2[N - 3][j] = X0[N - 2][j]

		opt2 = Cholesky(X2)
		precision2 = [1.0/s for s in error2]
		precision2[N-3] = 0.000001
		precision2[N-2] = precision0[N-2] + precision0[N-3]
		print(opt2, calc(X0, precision2))



		# if opt0 > (opt1 + opt2)/2.0:
		# if opt0 < max(opt1, opt2):
		if opt0 < min(opt1, opt2):
			print(opt0, opt1, opt2)
			print(X0)
			print(X1)
			print(X2)
			print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			break
