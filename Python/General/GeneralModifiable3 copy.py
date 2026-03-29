# General covariance
# The third and hopefully last version!

import numpy as np
import scipy
import random
import scipy.linalg
from pathlib import Path

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=200000, threshold=np.inf)

modifiableX = None
modifiableS = None
Z = None
P = 1.0
H = 0.5

FIXED_COORDS = False
FIXED_POINT = False



# ===================================================================

def get_v(N, X, S):
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

	return v


def get_X(v, N, X0, S0):
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

	return [N, X, S]


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

# The derivative of phi(r^2 / 2) = phi(x) with respect to x.
def dphi(x):
	return -np.exp(-x)

# Initial data is controlled via matrix Z and array modifiable.
def calc(args):
	global P, Z, H, modifiable, modifiableS, FIXED_POINT
	N = args[0]
	X = np.array(args[1])
	S = np.array(args[2])
	err2 = get_err2(N, S)

	D = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i + 1):
			D[i][j] = np.dot(X[i] - X[j], X[i] - X[j])
			D[j][i] = D[i][j]

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	dSigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(D[i][j] / 2.0)
			dSigma[i][j] = dphi(D[i][j] / 2.0)
		Sigma[i][i] += err2[i]

	C = scipy.linalg.cholesky(Sigma, lower=True)
	EVAL = np.dot(np.transpose(Z), X[N])

	if FIXED_POINT:
		ANS = np.dot(EVAL, C[N]) - 1.0 / (2.0 * H) * np.dot(X[N], X[N]) - H / 2.0
	else:
		ANS = np.dot(EVAL, C[N])
	print("ANS", ANS)
 

	# ======   Calculating gradient   ======

	K = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			K[i][j] = C[N][max(i, j)] * EVAL[min(i, j)]

	M = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(K), lower=False)
	M = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(M), lower=False)
	MdSigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			MdSigma[i][j] = M[i][j] * dSigma[i][j]

	ones_sums = MdSigma @ np.ones(N + 1)
	X_sums = MdSigma @ X


	gradient = []
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				partial = X[i][j] * ones_sums[i] -  X_sums[i][j]

				if i == N:
					partial += np.dot(Z[j], C[N])
					if FIXED_POINT:
						partial -= X[N][j] / H

				gradient.append(partial)

	for i in range(N):
		if modifiableS[i]:
			partial = MdSigma[i][i] * err2[i] / 2.0

			for j in range(N):
				if modifiableS[j]:
					partial -= MdSigma[j][j] * np.exp(S[i] - S[j]) / (2.0 * P)

			gradient.append(partial)


	ANS = -ANS
	gradient = -np.array(gradient)
	return (ANS, np.array(gradient))



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
			X[i][j] = 1.0/np.sqrt(j + 1)
			# X[i][j] = 1.0 - 2.0 * random.random()
			# X[i][j] = (1.0 - 2.0 * random.random())/5.0

	return X

def minimize():
	global FIXED_COORDS, modifiable, modifiableS, Z
	FIXED_COORDS = False
	X0 = None
	S0 = None
	N = None

	# Reading from input file
	if False:
		INPUT = read_from_file("in.txt")
		X0 = np.array(INPUT[0])
		N = len(X0) - 1

		modifiableS = [True for i in range(N)]
		S0 = get_S(N, INPUT[1])
	else:
		N = 200
		X0 = initial_guess(N)

		S0 = np.zeros((N), dtype='double')
		modifiableS = [True for i in range(N)]

	modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
	for i in range(N + 1):
		for j in range(i):
			modifiable[i][j] = True

	Z = np.eye(N + 1, dtype='double')
	
	fun = lambda v: calc(get_X(v, N, X0, S0))
	# fun(get_v(N, X0, S0))
	RES = scipy.optimize.minimize(fun, get_v(N, X0, S0), jac=True, options={"gtol": 0.0000000000000000000001}, method='BFGS')
	
	print("=================================================================")
	print(RES)

	(N, X1, S1) = get_X(RES.x, N, X0, S0)
	print(X1)
	print(get_err2(N, S1))

	with open("out.txt", "a") as f:
		f.write(str(X1))
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


