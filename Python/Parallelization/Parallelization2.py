# Covariance is exp(-||x - y||^2 / 2)
# This implementation checks certain hypotheses

import numpy as np
import scipy
import random

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

def get_v(N, D, X, modifiable):
	v = []

	for i in range(N + 1):
		for j in range(D):
			if modifiable[i][j]:
				v.append(X[i][j])

	return v


def get_X(N, D, v, X0, modifiable):
	X = np.zeros((N + 1, D), dtype='double')

	k = 0
	for i in range(N + 1):
		for j in range(D):
			if modifiable[i][j]:
				X[i][j] = v[k]
				k += 1
			else:
				X[i][j] = X0[i][j]

	return X


# ===============================================================================

# X is strictly lower diagonal
def calc(X, err2, PARALLEL):
	N = len(X) - 1
	D = len(X[0])


	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		Sigma[i][i] += err2[i]
	C = scipy.linalg.cholesky(Sigma, lower=True)

	v = scipy.linalg.solve_triangular(C[:N, :N], scipy.linalg.solve_triangular(C[:N, :N], Sigma[N, :N], lower=True), lower=True, trans=1)


	LATENCY = 0.0
	for i in range(N - PARALLEL, N):
	# for i in range(N - PARALLEL, N - PARALLEL + 1):
		W = np.zeros((D), dtype='double')
		for j in range(N):
			W += -(X[i] - X[j]) * Sigma[i][j] * v[j]
		W -= -(X[i] - X[N]) * Sigma[i][N]
		LATENCY -= np.dot(W, W) * np.exp(1.2 * np.dot(X[i] - X[N], X[i] - X[N]))


	ANS = 0.0
	for k in range(N - PARALLEL, N):
		ANS += C[N][k]**2

	return (ANS, v, LATENCY)


# =========================================================================================

def initial_guess(N, D):
	X = np.zeros((N + 1, D), dtype='double')
	for i in range(N + 1):
		for j in range(D):
			# X[i][j] = random.random()
			# X[i][j] = 0.5 + 0.5*random.random()
			X[i][j] = 2.5 * random.random()

	return X


def minimize():
	N = 4
	D = 2
	PARALLEL = 2

	X0 = initial_guess(N, D)

	err2 = np.zeros((N + 1), dtype='double')
	for i in range(N + 1):
		err2[i] = 2.0**(5.0 - 10.0*random.random())
	
	modifiable = [[False for j in range(D)] for i in range(N + 1)]
	for i in range(N - PARALLEL, N):
	# for i in range(N - PARALLEL, N - PARALLEL + 1):
		for j in range(D):
			modifiable[i][j] = True

	fun = lambda v: -calc(get_X(N, D, v, X0, modifiable), err2, PARALLEL)[0]

	RES = scipy.optimize.minimize(fun, get_v(N, D, X0, modifiable), options={"gtol": 0.0000000000000000000001})
	print(RES)
	X1 = get_X(N, D, RES.x, X0, modifiable)
	(ANS, v, LATENCY) = calc(X1, err2, PARALLEL)

	print(X1)
	print(ANS)
	print(v)
	print(LATENCY)

	# First hypothesis: that if all weights are positive, then all evaluation points must be equal
	if False:
		if LATENCY < -0.000001:
			return False

		POSITIVE = True
		for i in range(N - PARALLEL, N):
			if v[i] < 0.0001:
				return False
			
		DIFFERENT = True
		for i in range(N - PARALLEL, N):
			for j in range(i + 1, N):
				if np.dot(X1[i] - X1[j], X1[i] - X1[j]) < 0.01:
					DIFFERENT = False

		if not DIFFERENT:
			return False
		
		print("======================================================   HYPOTHESIS 1 FAIL!!!")
		print("ANS:\t", ANS)
		print("LATENCY:\t", LATENCY)
		print("v:\t", v)
		print("err2:\t", err2)
		print("X0 = ")
		print(X0)
		print("X1 = ")
		print(X1)
		return True

	# Second hypothesis: that if any weight is negative, then moving corresponding point to X[N] will increase value.
	# (As long as we have optimized w.r.t. that point and the weight is non-positive)
	if True:
		if v[N-PARALLEL] < -0.0001:
			X2 = get_X(N, D, RES.x, X0, modifiable)
			for j in range(D):
				X2[N-PARALLEL][j] = X1[N][j]
			(ANS2, v2, LATENCY2) = calc(X2, err2, PARALLEL)

			# print("v1:")
			# print(v)
			# print("v2:")
			# print(v2)
			# print("err2:")
			# print(err2)
			# print("X1=")
			# print(X1)
			# print("X2=")
			# print(X2)
			# print("diff =", ANS2 - ANS)
			if ANS2 < ANS:
				print("======================================================   HYPOTHESIS 2 FAIL!!!")
				print("ANS:", ANS, ANS2)
				print("LATENCY:", LATENCY, LATENCY2)
				print("v:", v, v2)
				print("err2:", err2)
				print("X0 = ")
				print(X0)
				print("X1 = ")
				print(X1)
				print("X2 = ")
				print(X2)
				return True

	return False


for i in range(1):
	print("iteration", i)
	if minimize():
		break



# HYPOTHESIS 1 POTENTIAL COUNTEREXAMPLE:
# err2:    [0.04594550     0.63146102      0.18434014      0.06685135      0.55993601      2.77421508      2.59155930      0.78186395      0.86442665     ]
# X1 = 
# [[1.47257355     0.97563975     ]
#  [1.70506970     1.42510606     ]
#  [1.34764652     1.59002184     ]
#  [0.85596951     1.85350814     ]
#  [0.85795130     0.72274638     ]
#  [2.12943190     0.44697025     ]
#  [-0.21741421    -0.05328052    ]
#  [-0.03729832    2.68529236     ]
#  [0.17679576     2.48284111     ]]
# PARALLEL = 3
