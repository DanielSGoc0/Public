# The covariance function is e^(-||x - y||^2 / 2)
# In this modification, the oracle returns f(x_{Alg}) + err * changer(f(x_{Alg}), ||Df(x_{Alg})||^2),
# where err is orthogonal and has fixed modulus and changer is some function.
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

def changer(x, y):
	# return np.exp((x-1)**2)
	# return (x-1.3)**2
	# return 1.0
	# return (4.0 - x**2)/4.0
	# return (x-2)**2/4.0
	return 1.0
	# return y
	# return (4.0 - x**2)/8.0 + y

# =========================================================================================

# A is lower triangular (not strictly)
def calc(A, stats):
	global err2
	# print("A = ", A)
	A = np.array(A, dtype='double')
	N = A.shape[0]

	# X describes the geometry of points Xn and is strictly lower triangular
	# X2 is a list of magnitudes squared of vectors in X.
	# F describes the means of F(Xn)/sqrt(dim).
	# DF describes the geometry of gradients dF(Xn) and is lower triangular
	# Sigma_ij = e^(-||Xi - Xj||^2 / 2)
	# DF2 is the modulus squared of the true derivative
	X = np.zeros((N + 1, N + 1), dtype='double')
	X2 = np.zeros((N + 1), dtype='double')
	F = np.zeros((N + 1), dtype='double')
	DF = np.zeros((N + 1, N + 1), dtype='double')
	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	DF2 = np.zeros((N + 1), dtype='double')

	# initial values
	DF2[0] = 1.0
	DF[0][0] = math.sqrt(1.0 + changer(F[0], DF2[0]) * err2)
	Sigma[0][0] = 1.0


	for n in range(1, N + 1):
		# At n-th step we know (X[0], ... , X[n - 1]), (DF[0], ... , DF[n - 1])
		# We wish to learn the points X[n] and DF[n].

		# Firstly, we update the geometry of X[n] and Sigma:
		# X[n] = A[n - 1][0] * DF[0] + ... + A[n - 1][n - 1] * DF[n - 1]
		for i in range(n):
			X[n] += A[n - 1][i] * DF[i]

		X2[n] = np.dot(X[n], X[n])

		for i in range(n + 1):
			Sigma[n][i] = math.exp(-(X2[n] + X2[i]) / 2.0 + np.dot(X[i], X[n]))
			Sigma[i][n] = Sigma[n][i]


		# Next, calculate the magnitude^2 of orthogonal vector.
		v = np.array(Sigma[n, :n])
		Sigma_partial = np.array(Sigma[:n, :n])
		for i in range(n):
			Sigma_partial[i][i] += err2*changer(F[i], DF2[i])
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
				M[a*n + i][a*n + i] += err2*changer(F[a], DF2[a])

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

		F[n] = np.dot(l, r)


		# Finally, update the vector DF[n] accordingly.
		for i in range(n):
			DF[n][i] = mean[i]
		DF2[n] = c + np.dot(DF[n], DF[n])
		DF[n][n] = math.sqrt(c + err2*changer(F[n], DF2[n]))


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

	fun = lambda v: -calc(get_A(N, v), True)

	RES = scipy.optimize.minimize(fun, get_v(N, A0), options={"gtol": 0.000000000000000000001})
	print(RES)
	calc(get_A(N, RES.x), True)


minimize(3)
# calc([[0.5, 0, 0], [0.3, 0.8, 0], [0.7, 0.2, 0.4]], True)

# F = 
# [0.00000000	 0.33786044	 0.56394258	 0.75120057	]
# A = 
# [[0.68985280	 0.00000000	 0.00000000	]
#  [0.71281127	 0.58573487	 0.00000000	]
#  [0.70826891	 0.57274675	 0.39860702	]]
# X = 
# [[0.00000000	 0.00000000	 0.00000000	 0.00000000	]
#  [1.19486009	 0.00000000	 0.00000000	 0.00000000	]
#  [1.16378981	 0.94122790	 0.00000000	 0.00000000	]
#  [1.15001227	 0.85477556	 0.60008926	 0.00000000	]]
# DF = 
# [[1.73205081	 0.00000000	 0.00000000	 0.00000000	]
#  [-0.12093445	 1.60691801	 0.00000000	 0.00000000	]
#  [-0.01876702	 -0.16452672	 1.50546587	 0.00000000	]
#  [0.00000001	 0.00000002	 0.00000000	 1.38090007	]]
# [[3.00000000	 -0.20946461	 -0.03250543	 0.00000002	]
#  [-0.20946461	 2.59681063	 -0.26211137	 0.00000002	]
#  [-0.03250543	 -0.26211137	 2.29384872	 0.00000000	]
#  [0.00000002	 0.00000002	 0.00000000	 1.90688499	]]

print(calc([[0, 0, 0, 0], [0.5, 0, 0, 0], [0.3, 0.8, 0, 0], [0.7, -0.2, 0.4, 0]], False))
