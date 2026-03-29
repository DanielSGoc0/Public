# The covariance function is e^(-||x - y||^2 / 2)
# The computation model is different.
# allows for applying error after evaluation.
# Also contains an option to change the formula!
import scipy
import numpy as np
import random

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

# =========================================================================================

def get_v(N, X):
	v = []
	for i in range(N + 1):
		for j in range(i):
			v.append(X[i][j])
	return v

def get_X(N, v):
	X = np.zeros((N + 1, N + 1), dtype='double')
	k = 0
	for i in range(N + 1):
		for j in range(i):
			X[i][j] = v[k]
			k += 1
	return X

# =========================================================================================

# X is strictly lower triangular
def calc(X):
	X = np.array(X, dtype='double')
	N = X.shape[0] - 1

	# Sigma_ij = e^(-||Xi - Xj||^2 / 2) + delta_{i, j} err2
	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)

	Chol = scipy.linalg.cholesky(Sigma, lower=True)

	RES = np.dot(X[N], Chol[N])
	# RES -= 1.0/4 * np.dot(X[N], X[N])
	# RES += Chol[N][N]**2
	# RES -= 1.0/3 * np.sqrt(np.dot(X[N], X[N]))
	# RES -= np.exp(np.dot(X[N], X[N]) / 2.0) / 2.0
	print(RES)
	return RES

# =========================================================================================

def random_matrix(N):
	X0 = np.zeros((N + 1, N + 1), dtype='double')

	for i in range(N + 1):
		for j in range(N + 1):
			if j < i:
				# X0[i][j] = 1 - 2*random.random()
				X0[i][j] = random.random()
			else: 
				X0[i][j] = 0.0

	return X0

def random_matrix2(N):
	X0 = np.zeros((N + 1, N + 1), dtype='double')
	for j in range(N):
		# X0[N][j] = 1 - 2*random.random()
		X0[N][j] = random.random()

	for i in range(N):
		for j in range(i):
			X0[i][j] = X0[N][j]

	return X0

def minimize(N):
	X0 = random_matrix(N)
	print(X0)

	fun = lambda v: -calc(get_X(N, v))

	RES = scipy.optimize.minimize(fun, get_v(N, X0), options = {'gtol': 0})
	X = get_X(N, RES.x)
	print(RES)
	print(X)

minimize(1)
# MAX = 0

# N = 4
# X_MAX = random_matrix(N)
# h = 0.001
# for k in range(1000000):
# 	X1 = random_matrix2(N)
# 	dx = np.array([1 - 2*random.random() for i in range(N)])
# 	dx /= np.sqrt(np.dot(dx, dx))
# 	X0 = np.zeros((N + 1, N + 1), dtype='double')
# 	X2 = np.zeros((N + 1, N + 1), dtype='double')
# 	for i in range(N + 1):
# 		for j in range(i):
# 			X0[i][j] = X1[i][j]
# 			X2[i][j] = X1[i][j]
# 	for j in range(N):
# 		X0[N][j] -= dx[j] * h
# 		X2[N][j] += dx[j] * h
# 	v0 = calc(X0)
# 	v1 = calc(X1)
# 	v2 = calc(X2)
# 	ddf = (v1 - (v0 + v2) / 2.0)/(h**2 * np.dot(dx, dx) / 2.0)

# 	print(ddf)
# 	if MAX < ddf:
# 		MAX = ddf
# 		X_MAX = X1
# 	print("MAX: ", MAX)
# 	print("X_MAX: ", X_MAX)


# [[0.00000000     0.00000000      0.00000000      0.00000000     ]
#  [0.55598160     0.00000000      0.00000000      0.00000000     ]
#  [0.84274604     -0.10405689     0.00000000      0.00000000     ]
#  [0.65164548     0.40564545      -0.25208198     0.00000000     ]]
