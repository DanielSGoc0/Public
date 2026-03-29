# Covariance is exp(-||x - y||^2 / 2)
# This implementation checks certain hypotheses
# This time we allow for changes in err2.

import numpy as np
import scipy
import random

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

def get_v(N, D, X, S, modifiable):
	v = []

	for i in range(N + 1):
		for j in range(D):
			if modifiable[i][j]:
				v.append(X[i][j])

	for s in S:
		v.append(s)

	return v


def get_X(N, D, v, X0, modifiable, PARALLEL):
	X = np.zeros((N + 1, D), dtype='double')

	k = 0
	for i in range(N + 1):
		for j in range(D):
			if modifiable[i][j]:
				X[i][j] = v[k]
				k += 1
			else:
				X[i][j] = X0[i][j]

	S = np.zeros((PARALLEL), dtype='double')
	for i in range(PARALLEL):
		S[i] = v[k]
		k += 1

	return (X, S)


# ===============================================================================

def calc(args, err2, TOTAL):
	X = np.array(args[0])
	S = args[1]
	N = len(X) - 1
	D = len(X[0])
	PARALLEL = len(S)

	SUM = 0.0
	for i in range(PARALLEL):
		SUM += np.exp(S[i])
	for i in range(PARALLEL):
		err2[N - PARALLEL + i] = 1.0/(np.exp(S[i]) * TOTAL / SUM) 


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

	# LATENCY -= ((v[N-2] * err2[N-2])**2 - (v[N-1] * err2[N-1])**2)**2

	Z = np.zeros((D, N - PARALLEL), dtype='double')
	Z[0][0] = 0.0
	# Z[0][1] = -30.0
	Z[0][1] = -0.2

	# This is the correct version
	F = 0.0
	for i in range(N-PARALLEL):
		for j in range(D):
			F += (X[N][j] - X[i][j]) * Z[j][i] * C[N][i]

	S = 0.0
	for k in range(N - PARALLEL, N):
		S += C[N][k]**2
	S = np.sqrt(S)

	r = (np.sqrt(F**2 + 4.0 * S**2) - F)/(2.0 * S)
	ANS = np.exp(-r**2 / 2.0) * (F + r * S)
	print(F, ANS - F)

	ANS -= 0.017750360869124942

	return (ANS, v, LATENCY)


# =========================================================================================

def initial_guess(N, D):
	X = np.zeros((N + 1, D), dtype='double')
	for i in range(N + 1):
		for j in range(D):
			# X[i][j] = random.random()
			# X[i][j] = 0.5 + 0.5*random.random()
			X[i][j] = 3.0 - 6 * random.random()

	return X


def minimize():
	N = 4
	D = 1
	PARALLEL = 2

	X0 = initial_guess(N, D)
	X0[0] = [-0.3]
	X0[1] = [0.3]
	# X0[2] = [-1.0]
	# X0[3] = [1.0]
	X0[4] = [0.0]

	S0 = np.zeros((PARALLEL), dtype='double')
	err2 = np.zeros((N + 1), dtype='double')
	for i in range(N):
		err2[i] = np.exp(3.0 - 6.0*random.random())
	err2 = [0.0, 0.0, 2.0, 2.0, 0.0]

	TOTAL = 0.0
	for i in range(PARALLEL):
		TOTAL += 1.0/err2[N-PARALLEL + i]
		# S0[i] = np.log(1.0/err2[i])
		S0[i] = 10.0*random.random()
	# TOTAL = 2.0
	TOTAL = 0.02


	modifiable = [[False for j in range(D)] for i in range(N + 1)]
	for i in range(N - PARALLEL, N):
	# for i in range(N - PARALLEL, N + 1):
	# # for i in range(N - PARALLEL, N - PARALLEL + 1):
		for j in range(D):
			modifiable[i][j] = True

	fun = lambda v: -calc(get_X(N, D, v, X0, modifiable, PARALLEL), err2, TOTAL)[0]

	RES = scipy.optimize.minimize(fun, get_v(N, D, X0, S0, modifiable), options={"gtol": 0.0000000000000000000001})
	(X1, S1) = get_X(N, D, RES.x, X0, modifiable, PARALLEL)
	(ANS, v, LATENCY) = calc((X1, S1), err2, TOTAL)

	print("X1:")
	print(X1)
	print("S1:")
	print(S1)
	print("err2:")
	print(err2)
	print("TOTAL:")
	print(TOTAL)
	print("ANS =", ANS)
	print("v =", v)
	print("LATENCY =", LATENCY)
	print("k =", v[N-1] * err2[N-1])

	# Zeroth hypothesis: that once we allow the changes in err2 and PARALLEL=2, then the only option is 
	# to have the evaluation points equal.
	if abs(S1[0] - S1[1]) < 10.0 and np.dot(X1[N-1] - X1[N-2], X1[N-1] - X1[N-2]) > 0.01 and LATENCY > -0.00000001:
		print("X1:")
		print(X1)
		print("S1:")
		print(S1)
		print("err2:")
		print(err2)
		print("TOTAL:")
		print(TOTAL)
		print("ANS =", ANS)
		print("v =", v)
		print("LATENCY =", LATENCY)
		return True


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
	if False:
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

# HYPOTHESIS 0 COUNTEREXAMPLE:
# X1:
# [[0.45264289     1.32353742     ]
#  [0.99556115     0.94652209     ]
#  [1.27436658     1.46556686     ]
#  [0.81716596     1.75982440     ]
#  [0.08465216     0.02912560     ]]
# err2:
# [0.19242321      0.16633315      0.15517412      0.21409772      9.55166604     ]
# ANS = 0.05816040874933458
# v = [0.29803519  0.41814558      -0.21703142     -0.15730041    ]
# LATENCY = -1.505905638507132e-14
