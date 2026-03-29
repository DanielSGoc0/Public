# Used to evaluate stuff with original covariance matrix.
# We allow different covariance functions.
import scipy
import numpy as np

# =========================================================================================

def kronecker(i, j):
	if i == j:
		return 1.0
	else:
		return 0.0

# The original covariance function
def phi0(x, y):
	# return np.exp(-np.dot(x - y, x - y) / 2.0)
	return 1.0/(1.0 + np.dot(x - y, x - y) / 2.0)

# First derivative (no sign changes!!) with respect to r^2/2
def phi1(x, y):
	# return -np.exp(-np.dot(x - y, x - y) / 2.0)
	return -1.0/(1.0 + np.dot(x - y, x - y) / 2.0)**2

# Second derivative with respect to r^2/2
def phi2(x, y):
	# return np.exp(-np.dot(x - y, x - y) / 2.0)
	return 2.0/(1.0 + np.dot(x - y, x - y) / 2.0)**3


# =========================================================================================

# X is lower triangular (not strictly)
def calc(X, err2):
	print("X =")
	print(X)
	X = np.array(X, dtype='double')
	N = X.shape[0]

	# X describes the geometry of points Xn and is strictly lower triangular
	# F describes the means of F(Xn)/sqrt(dim).
	# DF describes the geometry of gradients dF(Xn)/sqrt(dim) and is lower triangular
	# Sigma0, Sigma1, Sigma2 are filled with phi0, phi1 and phi2
	F = np.zeros((N), dtype='double')
	DF = np.zeros((N, N), dtype='double')
	# Sigma0 = np.array([[phi0(X[i], X[j]) for j in range(N)] for i in range(N)], dtype='double')
	Sigma1 = np.array([[-phi1(X[i], X[j]) + kronecker(i, j) * err2[i] for j in range(N)] for i in range(N)], dtype='double')
	Sigma2 = np.array([[phi2(X[i], X[j]) for j in range(N)] for i in range(N)], dtype='double')

	C = scipy.linalg.cholesky(Sigma1, lower=True)
	DF[0][0] = C[0][0]
	F[0] = 0.0

	for n in range(1, N):
		# At n-th step we know (X[0], ... , X[n - 1]), (DF[0], ... , DF[n - 1])
		# We wish to learn the points DF[n] and F[n].

		# We create the matrix and fill it with data.
		M = np.zeros((n*n, n*n), dtype='double')
		for a in range(n):
			for b in range(n):
				for i in range(n):
					for j in range(n):
						M[a*n + i][b*n + j] = (Sigma1[a][b] * kronecker(i, j) - Sigma2[a][b] * (X[a][i] - X[b][i]) * (X[a][j] - X[b][j]))

		d = np.zeros(n*n, dtype='double')
		for a in range(n):
			for i in range(n):
				d[a*n + i] = DF[a][i]
		r = scipy.linalg.solve(M, d, assume_a='pos')

		L = np.zeros((n, n*n), dtype='double')
		for b in range(n):
			for i in range(n):
				for j in range(n):
					L[i][b*n + j] = (Sigma1[n][b] * kronecker(i, j) - Sigma2[n][b] * (X[n][i] - X[b][i]) * (X[n][j] - X[b][j]))

		mean = np.dot(L, r)


		# As well as the mean of F(Xn).
		l = np.zeros(n*n, dtype='double')
		for b in range(n):
			for j in range(n):
				l[b*n + j] = Sigma1[n][b] * (X[n][j] - X[b][j])

		F[n] = np.dot(l, r)


		# Finally, update the vector DF[n] accordingly.
		for i in range(n):
			DF[n][i] = mean[i]
		DF[n][n] = C[n][n]

	print("DF =")
	print(DF)
	print("F = ")
	print(F)
	# return F[N]

# =========================================================================================

def calc2(X, err2):
	X = np.array(X, dtype='double')
	N = X.shape[0]

	Sigma = np.array([[-phi1(X[i], X[j]) + kronecker(i, j) * err2[i] for j in range(N)] for i in range(N)], dtype='double')
	C = scipy.linalg.cholesky(Sigma, lower=True)
	
	F = np.zeros((N), dtype='double')
	DF = np.zeros((N, N), dtype='double')

	# Firstly check DF
	for n in range(N):
		for i in range(n + 1):
			DF[n][i] = C[n][i]
			F[n] += X[n][i] * C[n][i]
	
	print("DF =")
	print(DF)
	print("F = ")
	print(F)





# =========================================================================================

X = [[0, 0, 0, 0], [0.5, 0, 0, 0], [0.3, 0.8, 0, 0], [0.7, -0.2, 0.4, 0]]
err2 = [1.0, 2.1, 4.0, 0.4]

calc(X, err2)
calc2(X, err2)
