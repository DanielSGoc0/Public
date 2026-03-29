# The covariance function is e^(-||x - y||^2 / 2)
# In this modification, the query points satisfy <xi - xj, xj> = 0 for all i > j.
import scipy
import numpy as np
import random

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=130)

# =========================================================================================

def kronecker(i, j):
	if i == j:
		return 1.0
	else:
		return 0.0


# =========================================================================================
err2 = 1.0


# A is lower triangular (not strictly)
def calc(A, stats):
	global err2

	print("A = ", A)
	A = np.array(A, dtype='double')
	N = A.shape[0]

	# X describes the geometry of points Xn and is strictly lower triangular
	# X2 is a list of magnitudes squared of vectors in X.
	# F describes the means of F(Xn)/sqrt(dim).
	# DF describes the geometry of gradients dF(Xn) and is lower triangular
	# Sigma_ij = e^(-||Xi - Xj||^2 / 2)
	X = np.zeros((N + 1, N + 1), dtype='double')
	X2 = np.zeros((N + 1), dtype='double')
	F = np.zeros((N + 1), dtype='double')
	DF = np.zeros((N + 1, N + 1), dtype='double')
	Sigma = np.zeros((N + 1, N + 1), dtype='double')

	# initial values
	DF[0][0] = np.sqrt(1.0 + err2)
	Sigma[0][0] = 1.0
	
	for i in range(N + 1):
		for j in range(N + 1):
			if j < i:
				X[i][j] = A[j]


	for n in range(1, N + 1):
		# At n-th step we know (X[0], ... , X[n - 1]), (DF[0], ... , DF[n - 1])
		# We wish to learn the points DF[n].

		X2[n] = np.dot(X[n], X[n])

		for i in range(n + 1):
			Sigma[n][i] = np.exp(-(X2[n] + X2[i]) / 2.0 + np.dot(X[i], X[n]))
			Sigma[i][n] = Sigma[n][i]


		# Next, calculate the magnitude^2 of orthogonal vector.
		v = np.array(Sigma[n, :n])
		Sigma_partial = np.array(Sigma[:n, :n])
		for i in range(n):
			Sigma_partial[i][i] += err2
		c = 1.0 + err2 - np.dot(v, scipy.linalg.solve(Sigma_partial, v, assume_a='pos'))


		# And then calculate the projection of DF[n] onto first n vectors.
		M = np.zeros((n*n, n*n), dtype='double')
		for a in range(n):
			for b in range(n):
				for i in range(n):
					for j in range(n):
						M[a*n + i][b*n + j] = Sigma[a][b] * (kronecker(i, j) - (X[a][i] - X[b][i]) * (X[a][j] - X[b][j]))
		for a in range(n):
			for i in range(n):
				M[a*n + i][a*n + i] += err2

		d = np.zeros(n*n, dtype='double')
		for a in range(n):
			for i in range(n):
				d[a*n + i] = DF[a][i]
		r = scipy.linalg.solve(M, d, assume_a='pos')
		print("r = ")
		print(r)
		# vec = np.zeros(n*n, dtype='double')
		# for a in range(n):
		# 	for i in range(n):
		# 		vec[a*n + i] = Sigma[n - 1][a] * (X[n - 1][i] - X[a][i])
		# print("vec A inverse:")
		# print(scipy.linalg.solve(M, vec, assume_a='pos'))

		L = np.zeros((n, n*n), dtype='double')
		for b in range(n):
			for i in range(n):
				for j in range(n):
					L[i][b*n + j] = Sigma[n][b] * (kronecker(i, j) - (X[n][i] - X[b][i]) * (X[n][j] - X[b][j]))

		mean = np.dot(L, r)


		# As well as the mean of F(Xn).
		l = np.zeros(n*n, dtype='double')
		for b in range(n):
			for j in range(n):
				l[b*n + j] = Sigma[n][b] * (X[n][j] - X[b][j])

		F[n] = np.dot(l, r)


		# Finally, update the vector DF[n] accordingly.
		for i in range(n):
			DF[n][i] = mean[i]
		DF[n][n] = np.sqrt(c)


	if stats:
		print("==================   RESULTS:   ==================")
		print("F = ")
		print(F)
		print("A = ")
		print(A)
		print("X = ")
		print(X)
		print("DF = ")
		print(DF)
		print("Sigma inverse = ")
		print(np.linalg.inv(Sigma))
		# for i in range(N + 1):
		# 	print(np.dot(DF[i], DF[i]))
		print("other...")
		for i in range(1, N):
			print(DF[i][i]**2, 1 - np.exp(-A[i - 1]**2))
		for i in range(2, N + 1):
			print(DF[i - 1][i - 2], - A[i - 1] * A[i - 2] * np.exp(-A[i - 2]**2)/np.sqrt(1 - np.exp(-A[i - 2]**2)))
		for i in range(2, N):
			print(A[i - 1] * (A[i - 1] * np.sqrt(1 - np.exp(-A[i - 2]**2)) - A[i]  * np.exp(-A[i - 1]**2)/np.sqrt(1 - np.exp(-A[i - 1]**2))))
		print("Formula for T:")
		for k in range(1, N):
			SUM = 0.0
			for i in range(k):
				SUM += Sigma[k][i] * DF[i][i] * A[i]
			print(DF[k][k - 1], np.exp(-A[k - 1]**2 / 2) * DF[k - 1][k - 1] - A[k - 1] * SUM)
		print("General formula:")
		for k in range(2, N):
			print(DF[k - 1][k - 1]/A[k - 1] + A[k] * np.exp(-A[k - 1]**2 / 2) / DF[k][k], A[k - 1] / DF[k - 1][k - 1] + np.exp(-A[k - 2]**2 / 2)/A[k - 2] * DF[k - 2][k - 2])
		print("other....")
		for n in range(N):
			print(F[n + 1], np.exp(-A[n]**2 / 2) * (DF[n][n] * A[n] + F[n]))

		print("DFs subtracted:")
		ei = np.zeros((N + 1))
		ei[1] = 1.0
		ei_1 = np.zeros((N + 1))
		ei_1[0] = 1.0
		print(DF[1] - np.exp(-A[0]**2 / 2.0) * DF[0] - DF[1][1] * ei + A[0] * F[1] * ei_1)
		for i in range(2, N + 1):
			ei = np.zeros((N + 1))
			ei[i] = 1.0
			ei_1 = np.zeros((N + 1))
			ei_1[i-1] = 1.0
			ei_2 = np.zeros((N + 1))
			ei_2[i-2] = 1.0
			print(DF[i] - np.exp(-A[i - 1]**2 / 2.0) * DF[i - 1] - DF[i][i] * ei + A[i - 1] * F[i] * ei_1 - 1/DF[i - 1][i - 1]*np.exp(-A[i-2]**2 - A[i-1]**2 / 2.0)*A[i - 1]*A[i - 2]*ei_2)


		# v = np.zeros(N + 1)
		# for i in range(N + 1):
		# 	v[i] = Sigma[i][N]
		# 	print("i = ", i)
		# 	print(np.linalg.solve(Sigma, v))

	# print(scipy.linalg.cholesky(M))

	return F[N]


# =========================================================================================



def minimize(N):
	A0 = np.zeros((N), dtype='double')
	for i in range(N):
		A0[i] = 1.0 - (i + 0.0)/N

	fun = lambda v: -calc(v, True)

	RES = scipy.optimize.minimize(fun, A0, options={"gtol": 0.00000000001})
	print("===================================================================")
	calc(RES.x, True)
	# print(RES)


minimize(8)

# A0 = [1.28418252,	 1.03621857,	 0.87020880,	 0.68875856,	 0.40949061]
# A0 = [1.28418252,	 1.23621857,	 0.47020880,	 0.98875856,	 0.60949061]
# A0 = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
# A0 = [0.0, 0.0, 0.0, 0.0, 0.0]
# calc(A0, True)
