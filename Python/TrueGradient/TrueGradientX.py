# The covariance function is e^(-||x - y||^2 / 2)
# No modifications. This problem has been essentially solved.
import scipy
import numpy as np
import random

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=130)

err2 = 10000.0

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


# =========================================================================================

# X is strictly lower triangular
def calc(X):
	global err2

	X = np.array(X, dtype='double')
	N = X.shape[0] - 1
	# print("X = ", X)

	# F describes the means of F(Xn)/sqrt(dim).
	# DF describes the geometry of gradients dF(Xn) and is lower triangular
	# Sigma_ij = e^(-||Xi - Xj||^2 / 2)
	F = np.zeros((N + 1), dtype='double')
	DF = np.zeros((N + 1, N + 1), dtype='double')
	Sigma = np.zeros((N + 1, N + 1), dtype='double')

	# initial values
	DF[0][0] = np.sqrt(1.0 + err2)
	Sigma[0][0] = 1.0


	for n in range(1, N + 1):
		# At n-th step we know (X[0], ... , X[n - 1]), (DF[0], ... , DF[n - 1])
		# We wish to learn the points X[n] and DF[n].

		# Firstly, we update the geometry of Sigma:
		for i in range(n + 1):
			Sigma[n][i] = np.exp(-np.dot(X[i] - X[n], X[i] - X[n]) / 2.0)
			Sigma[i][n] = Sigma[n][i]


		# Next, calculate the magnitude^2 of orthogonal vector.
		v = np.array(Sigma[n, :n])
		Sigma_partial = np.array(Sigma[:n, :n])
		for i in range(n):
			Sigma_partial[i][i] += err2
		c = 1 - np.dot(v, scipy.linalg.solve(Sigma_partial, v, assume_a='pos'))


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
		# print("r = ")
		# print(r)

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
				# l[b*n + j] = -Sigma[n][b] * X[b][j]

		# print(l)
		# print(r)

		F[n] = np.dot(l, r)


		# Finally, update the vector DF[n] accordingly.
		for i in range(n):
			DF[n][i] = mean[i]
		DF[n][n] = np.sqrt(c + err2)

	print(DF)
	print("F = ", F)
	# print("Geo = ", np.dot(DF, np.transpose(DF)))
	# G = [[np.dot(X[i] - X[j], X[i] - X[j]) for j in range(N + 1)] for i in range(N + 1)]
	# print(np.array(G))
	return F[N]


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

	fun = lambda v: -calc(get_A(N, v))

	RES = scipy.optimize.minimize(fun, get_v(N, A0))
	X = get_A(N, RES.x)
	print(RES)


# minimize(5)
# calc([[0, 0, 0, 0], [0.5, 0, 0, 0], [0.3, 0.8, 0, 0], [0.7, -0.2, 0.4, 0]])
calc([[0.00000000, 0.00000000, 0.00000000, 0.00000000], [1.30000000, 0.00000000, 0.00000000, 0.00000000], [0.60000000, 1.52464985, 0.00000000, 0.00000000], [1.47800272, -0.49821413, 0.82310486, 0.00000000]])
