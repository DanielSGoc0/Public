# The covariance function K(x, y) is func(||x - y||^2)
# This modification has different covariance function.
import scipy
import numpy as np
import math
import random

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=130)
err2 = 0.0

# =========================================================================================

def kronecker(i, j):
	if i == j:
		return 1.0
	else:
		return 0.0

def get_v(N, A):
	v = []
	for i in range(N):
		for j in range(i + 1):
			v.append(A[i][j])
	return v

def get_A(N, v):
	A = np.zeros((N, N), dtype='double')
	k = 0
	for i in range(N):
		for j in range(i + 1):
			A[i][j] = v[k]
			k += 1
	return A

# If ||x - y||^2 = m2, then returns K(x, y)
def func(m2):
	# return np.exp(-m2/2)
	# return np.exp(-m2)
	return np.exp(-m2/2) + 0.1 * np.exp(-20*m2/2)

# If ||x - y||^2 = m2 and <x - y, ei> = a then returns D K(x, y) (ei)
def dfunc(m2, a, i):
	# return -np.exp(-m2/2) * a
	# return - 2 * np.exp(-m2) * a
	return -np.exp(-m2/2) * a - 0.1 * 20 * np.exp(-20*m2/2) * a

# If ||x - y||^2 = m2 and <x - y, ei> = a and <x - y, ej> = b then returns D^2 K(x, y) (ei, ej)
def ddfunc(m2, a, b, i, j):
	# return np.exp(-m2/2) * (kronecker(i, j) - a * b)
	# return np.exp(-m2) * (2*kronecker(i, j) - 4*a*b)
	return np.exp(-m2/2) * (kronecker(i, j) - a * b) + 0.1 * np.exp(-20*m2/2) * (20*kronecker(i, j) - 20**2*a*b)

# =========================================================================================

# A is lower triangular (not strictly)
def calc(A, stats):
	global err2
	print("A = ", A)
	A = np.array(A, dtype='double')
	N = A.shape[0]

	# X describes the geometry of points Xn and is strictly lower triangular
	# F describes the means of F(Xn)/sqrt(dim).
	# DF describes the geometry of gradients dF(Xn) and is lower triangular
	# Sigma_ij describes the diagonal entries of derivatives covariance matrices.
	X = np.zeros((N + 1, N + 1), dtype='double')
	F = np.zeros((N + 1), dtype='double')
	DF = np.zeros((N + 1, N + 1), dtype='double')
	Sigma = np.zeros((N + 1, N + 1), dtype='double')

	# initial values
	d_var = ddfunc(0, 0, 0, 0, 0)
	DF[0][0] = math.sqrt(d_var + err2)
	Sigma[0][0] = d_var

	for n in range(1, N + 1):
		# At n-th step we know (X[0], ... , X[n - 1]), (DF[0], ... , DF[n - 1])
		# We wish to learn the points X[n] and DF[n].

		# Firstly, we update the geometry of X[n] and Sigma:
		# X[n] = A[n - 1][0] * DF[0] + ... + A[n - 1][n - 1] * DF[n - 1]
		for i in range(n):
			X[n] += A[n - 1][i] * DF[i]

		for i in range(n + 1):
			Sigma[n][i] = ddfunc(np.dot(X[n] - X[i], X[n] - X[i]), 0, 0, 0, 0)
			Sigma[i][n] = Sigma[n][i]


		# Next, calculate the magnitude^2 of orthogonal vector.
		v = np.array(Sigma[n, :n])
		Sigma_partial = np.array(Sigma[:n, :n])
		for i in range(n):
			Sigma_partial[i][i] += err2
		c = d_var + err2 - np.dot(v, scipy.linalg.solve(Sigma_partial, v, assume_a='pos'))


		# And then calculate the projection of DF[n] onto first n vectors.
		M = np.zeros((n*n, n*n), dtype='double')
		for a in range(n):
			for b in range(n):
				m2 = np.dot(X[a] - X[b], X[a] - X[b])
				for i in range(n):
					for j in range(n):
						M[a*n + i][b*n + j] = ddfunc(m2, X[a][i] - X[b][i], X[a][j] - X[b][j], i, j)
		for a in range(n):
			for i in range(n):
				M[a*n + i][a*n + i] += err2

		d = np.zeros(n*n, dtype='double')
		for a in range(n):
			for i in range(n):
				d[a*n + i] = DF[a][i]
		r = scipy.linalg.solve(M, d, assume_a='pos')

		L = np.zeros((n, n*n), dtype='double')
		for b in range(n):
			m2 = np.dot(X[n] - X[b], X[n] - X[b])
			for i in range(n):
				for j in range(n):
					L[i][b*n + j] = ddfunc(m2, X[n][i] - X[b][i], X[n][j] - X[b][j], i, j)

		mean = np.dot(L, r)


		# As well as the mean of F(Xn).
		l = np.zeros(n*n, dtype='double')
		for b in range(n):
			m2 = np.dot(X[n] - X[b], X[n] - X[b])
			for j in range(n):
				l[b*n + j] = dfunc(m2, X[b][j] - X[n][j], j)

		F[n] = np.dot(l, r)


		# Finally, update the vector DF[n] accordingly.
		for i in range(n):
			DF[n][i] = mean[i]
		DF[n][n] = math.sqrt(c)


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
		print(np.dot(DF, np.transpose(DF)))

	return F[N]
	# return -np.dot(DF[N], DF[N])


# =========================================================================================

def initial_guess(N):
	A = np.zeros((N, N), dtype='double')

	for i in range(N):
		for j in range(N):
			if j <= i:
				A[i][j] = random.random()
			else: 
				A[i][j] = 0.0

	return A


def minimize(N):
	A0 = initial_guess(N)
	print(A0)

	fun = lambda v: -calc(get_A(N, v), False)

	RES = scipy.optimize.minimize(fun, get_v(N, A0), options={"gtol": 0.000000000000000000001})
	calc(get_A(N, RES.x), True)


minimize(3)
