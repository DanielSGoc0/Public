# Program created to analyze the parallelization hypothesis

import numpy as np
import random 

import scipy.linalg

N = 3
D = 1
error2 = np.zeros((N), dtype='double')

def Cholesky(X):
	global N, error2

	Sigma = np.zeros((N, N), dtype='double')
	for i in range(N):
		for j in range(N):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		Sigma[i][i] += error2[i]
		# Sigma[i][i] += 1.0

	C = np.linalg.cholesky(Sigma)
	print("===============")
	print(X)
	print(Sigma)
	print(error2)
	print(C[N - 1][N - 3]**2 + C[N - 1][N - 2]**2)
	v12 = Sigma[N-1][N-1] - np.dot(Sigma[N-1][:(N-1)], scipy.linalg.solve(Sigma[:(N-1), :(N-1)], Sigma[N-1][:(N-1)], assume_a='pos'))
	v1 = Sigma[N-1][N-1] - np.dot(Sigma[N-1][:(N-3)], scipy.linalg.solve(Sigma[:(N-3), :(N-3)], Sigma[N-1][:(N-3)], assume_a='pos'))
	print(v1 - v12)
	print(v1, v12)

	# return np.linalg.cholesky(Sigma)
	return v12

def optimal(X):
	global N, error2

	# C = Cholesky(X)
	# print(C)

	print("optimal N:", N)
	# print(C[N - 1][N - 3]**2 + C[N - 1][N - 2]**2)
	# return C[N - 1][N - 3]**2 + C[N - 1][N - 2]**2
	return Cholesky(X)

for s in range(100000):
# for s in range(1):
	if s % 1000 == 0:
		print("s =", s)
	error2 = [random.random() for k in range(N)]
	# error2 = [0.1 for k in range(N)]
	X0 = np.array([[random.random()*2.0 for j in range(D)] for i in range(N)])
	X0 = [[1.49161777], [0.37608213], [1.08462918]]
	X0 = np.array(X0)
	# X0[0][0] = 0.0
	# X0[0][1] = 0.0
	# X0[1][1] = 0.0

	opt0 = optimal(X0)


	# Chol = Cholesky(X0)[(N - 3):N, (N - 3):N]
	# M = np.dot(Chol, np.transpose(Chol))
	# for i in range(3):
	# 	M[i][i] -= error2[i + N - 3]
	# 	for j in range(3):
	# 		M[i][j] *= np.exp(np.dot(X0[i + N - 3] - X0[j + N - 3], X0[i + N - 3] - X0[j + N - 3]) / 2.0)
	# 		if M[i][j] > 1:
	# 			print(X0)
	# 			print(M)

	# print(M)



	if True:
		X1 = np.zeros((N, D), dtype='double')
		for i in range(N):
			for j in range(D):
				X1[i][j] = X0[i][j]
		for j in range(D):
			X1[N - 2][j] = X0[N - 3][j]

		opt1 = optimal(X1)

		X2 = np.zeros((N, D), dtype='double')
		for i in range(N):
			for j in range(D):
				X2[i][j] = X0[i][j]
		for j in range(D):
			X2[N - 3][j] = X0[N - 2][j]

		opt2 = optimal(X2)

		# X3 = np.zeros((N, D), dtype='double')
		# for i in range(N):
		# 	for j in range(D):
		# 		X3[i][j] = X0[i][j]
		# for j in range(D):
		# 	X3[N - 2][j] = X0[N - 1][j]
		# 	X3[N - 3][j] = X0[N - 1][j]

		# opt3 = optimal(X3)

		# X4 = np.zeros((N, D), dtype='double')
		# for i in range(N):
		# 	for j in range(D):
		# 		X4[i][j] = X0[i][j]
		# for j in range(D):
		# 	X4[N - 2][j] = (X0[N - 3][j] + X0[N - 2][j])/2.0
		# 	X4[N - 3][j] = (X0[N - 3][j] + X0[N - 2][j])/2.0
		# 	# X4[N - 2][j] = (1*X0[N - 3][j] + 2*X0[N - 2][j])/3.0
		# 	# X4[N - 3][j] = (1*X0[N - 3][j] + 2*X0[N - 2][j])/3.0

		# opt4 = optimal(X4)

		# X5 = np.zeros((N, D), dtype='double')
		# for i in range(N):
		# 	for j in range(D):
		# 		X5[i][j] = X0[i][j]
		# for j in range(D):
		# 	# X5[N - 2][j] = (X0[N - 3][j] + X0[N - 2][j])/2.0
		# 	# X5[N - 3][j] = (X0[N - 3][j] + X0[N - 2][j])/2.0
		# 	X5[N - 2][j] = (2*X0[N - 3][j] + 1*X0[N - 2][j])/3.0
		# 	X5[N - 3][j] = (2*X0[N - 3][j] + 1*X0[N - 2][j])/3.0

		# opt5 = optimal(X5)


		# if opt0 < (opt1 + opt2)/2.0:
		# if opt0 > (opt1 + opt2)/2.0:
		if opt0 > max(opt1, opt2):
			print(opt0, opt1, opt2)
			print(error2)
			print(X0)
			print(X1)
			print(X2)
			print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			break
		# if opt0 > (opt3 + opt4)/2.0:
		# 	print(opt0, opt3, opt4)
		# 	print("...")
		# if opt0 > max(opt3, opt4):
		# 	print(opt0, opt3, opt4)
		# 	print(X0)
		# 	print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		# if opt0 > max(opt3, opt4, opt5):
		# 	print(opt0, opt3, opt4, opt5)
		# 	print(X0)
		# 	print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		# if opt0 > max(opt1, opt2, opt3, opt4):
		# 	print(opt0, opt1, opt2, opt3, opt4)
		# 	print(X0)
		# 	print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
