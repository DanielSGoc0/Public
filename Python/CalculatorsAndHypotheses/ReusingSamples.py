# We check whether using samples (i.e. error coordinates)
# can be helpful or not. 

# It turns out that yes, it can be helpful
# However, that requires very strong assumptions, namely
# that we know the covariances between stuff, and that's
# too restrictive!

import numpy as np
import scipy
import scipy.linalg
from pathlib import Path
import random

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f} \t'.format}, linewidth=200000, threshold=np.inf)

modifiableX = None
modifiableS = None
P = 1.0
H = 0.5
R2 = None

FIXED_COORDS = False
FIXED_POINT = False



# ===================================================================

def get_v(N, M, X, W, S):
	global FIXED_COORDS, modifiable, modifiableS

	v = []
	if FIXED_COORDS:
		v = X[N][:N]
	else:
		for i in range(N + 1):
			for j in range(N + 1):
				if modifiable[i][j]:
					v.append(X[i][j])

	for i in range(len(modifiableS)):
		if modifiableS[i]:
			v.append(S[i])

	for i in range(M, N):
		for j in range(M, i):
			v.append(W[i][j])

	return v


def get_args(v, N, M, X0, S0):
	global FIXED_COORDS, modifiable, modifiableS

	X = np.zeros((N + 1, N + 1), dtype='double')
	k = 0
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				if FIXED_COORDS:
					X[i][j] = v[j]
				else:
					X[i][j] = v[k]
					k += 1
			else:
				X[i][j] = X0[i][j]
	if FIXED_COORDS:
		k = N

	S = np.zeros(len(modifiableS), dtype='double')
	for i in range(len(modifiableS)):
		if modifiableS[i]:
			S[i] = v[k]
			k += 1
		else:
			S[i] = S0[i]

	W = np.zeros((N, N), dtype='double')

	for i in range(M, N):
		for j in range(M, i):
			W[i][j] = v[k]
			k += 1
		W[i][i] = 1.0

	return [N, M, X, W, S]


def get_err2(N, S):
	global P, modifiableS

	err2 = np.zeros((N + 1), dtype='double')
	SUM = 0.0
	for i in range(N):
		if modifiableS[i]:
			SUM += np.exp(S[i])
	for i in range(N):
		if modifiableS[i]:
			err2[i] = 1.0/(np.exp(S[i]) * P / SUM)
		else:
			err2[i] = S[i]

	return err2


def get_S(N, err2):
	global P, modifiableS

	S = np.zeros((N), dtype='double')
	for i in range(N):
		if modifiableS[i]:
			S[i] = -np.log(P * err2[i])
		else:
			S[i] = err2[i]

	return S


# ===============================================================================

# The original covariance function phi(r^2 / 2).
# By Bernstein's Theorem, it must be an integral of exponentials.
def phi(x):
	return np.exp(-x)

def calc(args):
	global P, H, modifiable, modifiableS, FIXED_POINT
	N = args[0]
	OFFSET = args[1]
	X = np.array(args[2])
	W = np.array(args[3])
	S = np.array(args[4])
	err2 = get_err2(N, S)

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)

	# W = np.eye(N, dtype='double')
	# G = W @ W.T
	# f = np.zeros((N), dtype='double')
	# for i in range(OFFSET, N):
	# 	f[i] = np.sqrt(err2[i] / G[i][i])
	# for i in range(OFFSET, N):
	# 	for j in range(OFFSET, N):
	# 		G[i][j] *= f[i] * f[j]
	# 		Sigma[i][j] *= (1.0 + G[i][j])
	# print(G)

	# for i in range(OFFSET, N):
	# 	sum = 0.0
	# 	for j in range(i, OFFSET - 1, -1):
	# 		Sigma[i][j] += np.exp(-sum / 2.0)
	# 		if i != j:
	# 			Sigma[j][i] += np.exp(-sum / 2.0)
	# 		sum += 1.0 / err2[j]

	# for i in range(OFFSET, N, 2):
	# 	for j in range(OFFSET, N, 2):
	# 		Sigma[i][j] *= (1.0 + np.sqrt(err2[0] * err2[0]))
	# 		Sigma[i + 1][j + 1] *= (1.0 + np.sqrt(err2[1] * err2[1]))
	# 		# Sigma[i][j] *= (1.0 + np.sqrt(err2[i] * err2[j]))
	# 		# Sigma[i + 1][j + 1] *= (1.0 + np.sqrt(err2[i + 1] * err2[j + 1]))

	for i in range(OFFSET, N):
		Sigma[i][i] += err2[i % 2]

	C = scipy.linalg.cholesky(Sigma, lower=True)

	print(Sigma)
	print(C)

	ANS = np.dot(X[N, OFFSET:N], C[N, OFFSET:N])
	print(ANS)
	return -ANS


# =========================================================================================

def read_number(txt, cursor):
	x = 0
	negative = False
	after_decimal = -1
	while True:
		c = txt[cursor]
		if c == ' ' or c == '\t' or c == '\n' or c == ']':
			break
		elif c == '-':
			cursor += 1
			negative = True
		elif c == '.':
			cursor += 1
			after_decimal = 0
			continue
		else:
			cursor += 1
			x *= 10
			x += ord(c) - ord('0')
		if after_decimal >= 0:
			after_decimal += 1
	if after_decimal >= 0:
		x = float(x)
		x /= 10**after_decimal
	if negative:
		x = -x
	return (x, cursor)

def create_array(txt, cursor):
	RES = []
	while True:
		if cursor == len(txt):
			return (RES, cursor)
		c = txt[cursor]
		if c == ' ' or c == '\t' or c == '\n':
			cursor += 1
			continue
		elif c == '[':
			cursor += 1
			output = create_array(txt, cursor)
			RES.append(output[0])
			cursor = output[1]
		elif c == ']':
			cursor += 1
			return (RES, cursor)
		else:
			output = read_number(txt, cursor)
			RES.append(output[0])
			cursor = output[1]

def read_from_file(filename):
	txt = Path(filename).read_text()
	RES = create_array(txt, 0)[0]
	return RES


# ============================================================================


def initial_guess(N):
	X = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i):
			# X[i][j] = 1.0/np.sqrt(j + 1)
			# X[i][j] = 0.01
			# X[i][j] = 1.0/np.sqrt(j + 1) + random.random() / 10.0
			# X[i][j] = 0.071
			# X[i][j] = 1.0 - random.random()
			X[i][j] = (1.0 - 2.0 * random.random())
			# if j == 0:
			# 	X[i][j] = (-1)**i
			# elif j == 1:
			# 	X[i][j] = 0.0

	return X

def minimize():
	global FIXED_COORDS, modifiable, modifiableS, R2, P
	FIXED_COORDS = False
	X0 = None
	W0 = None
	S0 = None
	N = None
	M = None

	# Reading from input file
	# Indices 0 to M - 1 are meant for initial data.
	# Indices M to N - 1 are meant for evaluation points.
	# Index N is meant for output.
	if False:
		INPUT = read_from_file("in.txt")
		X0 = np.array(INPUT[0])
		N = len(X0) - 1
		Z0 = np.array(INPUT[2])
		M = len(Z0)

		modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
		for i in range(M, N + 1):
			for j in range(i):
				modifiable[i][j] = True

		modifiableS = [True for i in range(N)]
		for i in range(M):
			modifiableS[i] = False
		for i in range(M, N):
			modifiableS[i] = False
		S0 = get_S(N, INPUT[1])

		Z = np.eye(N + 1, dtype='double')
		for i in range(M):
			for j in range(M):
				Z[i][j] = Z0[i][j]
	else:
		N = 4
		M = 0
		P = 0.5

		X0 = initial_guess(N)
		# X0[0][0] = 0.5798527
		# X0[0][1] = 0.32425464266
		# X0[0][2] = 0.163347653754
		# X0[1][0] = -0.57285235
		# X0[1][1] = -0.067735435654
		# X0[1][2] = -0.35463773
		# X0[2][0] = 0.867143785681
		# X0[2][1] = 0.538572989
		# X0[2][2] = -0.14325635
		# X0[N][0] = 0.0
		# X0[N][1] = 0.0
		# X0[N][2] = 0.0

		# X0[0][0] = -0.3
		# X0[0][1] = 0.0
		# X0[1][0] = 0.3
		# X0[1][1] = 0.0
		# X0[N][0] = 0.0
		# X0[N][1] = 0.0

		S0 = np.ones((N), dtype='double')
		# S0[0] = 4.418927941
		# S0[1] = 0.834729724
		# S0[2] = 2.943724987

		# S0[0] = 0.0
		# S0[1] = 0.0

		modifiableS = [True for i in range(N)]
		for i in range(M):
			modifiableS[i] = False
		for i in range(M, N):
			# modifiableS[i] = False
			S0[i] = (N - M)/P
	
		# Z[1][0] = -0.58743289793
		# Z[2][0] = 0.252353
		# Z[1][2] = 1.343151533
		# Z[2][1] = -0.375734553
		# Z[0][1] = 0.83925782935
		# Z[0][2] = -0.2342352359

		# Z[0][0] = 0.0
		# Z[0][1] = -2.0


		modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
		for i in range(M, N + 1):
			for j in range(i):
				modifiable[i][j] = True
		for i in range(M + 1, N, 2):
			modifiable[i][i - 1] = False
			X0[i][i - 1] = 0.0
		for j in range(M):
			modifiable[N][j] = False
		# for i in range(N + 1):
		# 	modifiable[i][1] = False
		# 	X0[i][1] = 0.0

		# for i in range(3, N, 2):
		# 	modifiable[i][i - 1] = False
		# 	X0[i][i - 1] = 0.0
		
		# R2 = 6.8967142387450595
		# X0 = initial_state(N, M, X0)

		W0 = np.eye(N, dtype='double')

	# print(X0)
	# print()
	# err0 = get_err2(N, S0)
	# while True:
	# # for i in range(2):
	# 	(X0, err0) = iterate(N, M, X0, err0)
	# 	time.sleep(0.1)
	# 	# print()
	# return
		
	
	fun = lambda v: calc(get_args(v, N, M, X0, S0))

	# (ANS0, grad) = fun(get_v(N, X0, S0))
	# for i in range(N - 1):
	# 	X1 = initial_guess(N)
	# 	X2 = initial_guess(N)
	# 	S1 = np.copy(S0)
	# 	S2 = np.copy(S0)
	# 	# X1[N - 1][i] += 0.001
	# 	# X2[N - 1][i] -= 0.001
	# 	X1[N][i] += 0.001
	# 	X2[N][i] -= 0.001
	# 	# S1[M + i] += 0.001
	# 	# S2[M + i] -= 0.001

	# 	(ANS1, _1) = fun(get_v(N, X1, S1))
	# 	(ANS2, _2) = fun(get_v(N, X2, S2))
	# 	# print((ANS1 - ANS2) / 0.002, grad[((N - 2) * (N - 1) // 2) - (M * (M - 1) // 2) + i])
	# 	print((ANS1 - ANS2) / 0.002, grad[(N * (N - 1) // 2) - (M * (M - 1) // 2) + i])
	# 	# print((ANS1 - ANS2) / 0.002, grad[(N * (N + 1) // 2) - (M * (M - 1) // 2) + i])


	# fun(get_v(N, X0, S0))
	# for j in range(M):
	# 	modifiable[N][j] = False
	RES = scipy.optimize.minimize(fun, get_v(N, M, X0, W0, S0), options={"gtol": 0.0000000000000000000001}, method='BFGS')
	
	print("=================================================================")
	print(RES)

	(N, M, X1, W1, S1) = get_args(RES.x, N, M, X0, S0)
	print(X1)
	print(W1)
	print(get_err2(N, S1))

	with open("out.txt", "a") as f:
		f.write(str(RES))
		f.write("\n")
		f.write(str(X1))
		f.write("\n")
		f.write(str(W1))
		f.write("\n")
		f.write(str(get_err2(N, S1)))
		f.write("\n\n\n")

	# return (fun(RES.x)[0], X1)


minimize()

# best_ans = 100000.0
# best_X = None

# for k in range(10):
# 	print(k)
# 	(ans, X) = minimize()
# 	print(ans)

# 	if ans < best_ans:
# 		best_ans = ans
# 		best_X = X

# print(best_ans)
# print(best_X)
