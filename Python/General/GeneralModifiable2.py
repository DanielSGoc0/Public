# Covariance is exp(-||x - y||^2 / 2)
# A most general implementation of Maximization problem/Fixed Point Problem.
# The second version, with much more efficient gradient!

import numpy as np
import scipy
import random
import scipy.linalg
from pathlib import Path

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=200000, threshold=np.inf)

modifiable = None
modifiableS = None
Z = None
P = 1.0
FIXED_COORDS = False

H = 0.5
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


def get_X(N, v, X0, S0):
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

	return [X, S]


def get_err2(S):
	global P, modifiableS
	N = len(modifiableS)

	err2 = np.zeros((N + 1), dtype='double')
	SUM = 0.0
	for i in range(N):
		if modifiableS[i]:
			SUM += np.exp(S[i])
	for i in range(N):
		if modifiableS[i]:
			err2[i] = 1.0/(np.exp(S[i]) * P / SUM)

	return err2


def get_S(err2):
	global P, modifiableS
	N = len(err2) - 1

	S = np.zeros((N), dtype='double')
	for i in range(N):
		if modifiableS[i]:
			S[i] = -np.log(P * err2[i])

	return S


# ===============================================================================

# X is strictly lower diagonal
def calc(args):
	global P, Z, modifiable, modifiableS, H, FIXED_POINT
	X = args[0]
	S = args[1]
	N = len(X) - 1
	X = np.array(X)

	err2 = np.zeros((N + 1), dtype='double')
	SUM = 0.0
	for i in range(N):
		if modifiableS[i]:
			SUM += np.exp(S[i])
	for i in range(N):
		if modifiableS[i]:
			err2[i] = 1.0/(np.exp(S[i]) * P / SUM)
	# err2 = np.ones((N + 1), dtype='double')
	# for i in range(N):
	# 	err2[i] = 10.0
	# err2[N] = 100000.0

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
		Sigma[i][i] += err2[i]

	C = scipy.linalg.cholesky(Sigma, lower=True)
	CinvT = np.transpose(scipy.linalg.solve_triangular(C, np.eye(N + 1), lower=True))

	EVAL = np.dot(np.transpose(Z), X[N])

	if FIXED_POINT:
		ANS = np.dot(EVAL, C[N]) - 1.0/(2.0 * H) * np.dot(X[N], X[N]) - H / 2.0
	else:
		ANS = np.dot(EVAL, C[N])
 



	Isigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N):
		Isigma[i][i] = np.sqrt(err2[i])
	Xsigma = np.dot(X, Isigma)
	Csigma = np.dot(C, Isigma)

	K = np.zeros((N + 1, N + 1), dtype='double')
	Kprime = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			K[i][j] = C[N][max(i, j)] * X[N][min(i, j)]
			Kprime[i][j] = X[N][max(i, j)] * C[N][min(i, j)]

	M = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(K), lower=False)
	M = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(M), lower=False)
	MSigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			MSigma[i][j] = M[i][j] * Sigma[i][j]

	XNXN = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			XNXN[i][j] = X[N][i] * X[N][j]

	# print(M)
	# print(MSigma)
	# print(Xsigma)
	# print(np.dot(np.transpose(Xsigma), np.dot(MSigma, Xsigma)))







	gradient = []
	# correlations = []

	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				# Append gradient wrt single variable Xij
				partial = 0.0

				v = np.zeros((N + 1), dtype='double')
				for s in range(N + 1):
					v[s] = -(X[i][j] - X[s][j]) * Sigma[i][s]

				a = scipy.linalg.solve_triangular(C, v, lower=True)
				b = CinvT[i]
				# w = np.zeros((N + 1), dtype='double')

				# Apply phi to a b^T + b a^T
				summand = 0.0
				for s in range(N + 1):
					partial += summand * C[N][s] * a[s]
					summand += EVAL[s] * b[s]

				# summand = 0.0
				# for s in range(N, -1, -1):
				# 	w[s] += b[s] * summand
				# 	summand += C[N][s] * a[s]

				summand = 0.0
				for s in range(N + 1):
					partial += summand * C[N][s] * b[s]
					summand += EVAL[s] * a[s]

				# summand = 0.0
				# for s in range(N, -1, -1):
				# 	w[s] += a[s] * summand
				# 	summand += C[N][s] * b[s]

				for s in range(N + 1):
					partial += C[N][s] * EVAL[s] * a[s] * b[s]
					
				# for s in range(N + 1):
				# 	w[s] += C[N][s] * a[s] * b[s]

				# We've calculated the ingredient coming from dCN / dxij. 
				# If i = N, we also have to append the part coming from dxN
				if i == N:
					partial += np.dot(Z[j], C[N])

					if FIXED_POINT:
						partial -= X[N][j] / H

				gradient.append(partial)
				# if i != N:
				# 	# correlations.append(w[N])
				# 	w[N] = 0.0
				# 	# correlations.append(np.log(np.dot(w, w)))
				# 	correlations.append(np.dot(w, X[N]) / np.sqrt(np.dot(X[N], X[N]) * np.dot(w, w)))
					# correlations.append((partial - np.dot(w, X[N]))/(np.abs(partial) + np.abs(np.dot(w, X[N]))))

	if len(modifiableS) > 0:
		# We also have to deal with change of errors.
		# Firstly calculate the gradient wrt change of a single diagonal element
		diagonals = np.zeros((N), dtype='double')
		# W = np.zeros((N, N + 1), dtype='double')

		for i in range(N):
			partial = 0.0
			b = CinvT[i]

			# Apply phi to b b^T
			summand = 0.0
			for s in range(N + 1):
				partial += summand * C[N][s] * b[s]
				summand += EVAL[s] * b[s]
			
			# summand = 0.0
			# for s in range(N, -1, -1):
			# 	W[i][s] += b[s] * summand
			# 	summand += C[N][s] * b[s]

			for s in range(N + 1):
				partial += C[N][s] * EVAL[s] * b[s] * b[s] / 2.0
			
			# for s in range(N + 1):
			# 	W[i][s] += C[N][s] * b[s] * b[s] / 2

			diagonals[i] = partial


		for i in range(len(modifiableS)):
			if modifiableS[i]:
				# Append gradient wrt single variable Si.

				partial = -diagonals[i] * err2[i]
				# w = -W[i] * err2[i]

				for j in range(N):
					if modifiableS[j]:
						partial += diagonals[j] * np.exp(S[i] - S[j]) / P
				
				# for j in range(N):
				# 	if modifiableS[j]:
				# 		w += W[j] * np.exp(S[i] - S[j]) / P

				gradient.append(partial)
				# correlations.append(np.log(np.dot(w, w)))
				# correlations.append(np.dot(w, X[N]) / np.sqrt(np.dot(X[N], X[N]) * np.dot(w, w)))
				# correlations.append((partial - np.dot(w, X[N]))/(np.abs(partial) + np.abs(np.dot(w, X[N]))))

	print(gradient)
	print("ANS", ANS)
	# print(np.array(correlations))

	# print(C)
	# print(np.linalg.solve(C, np.eye(N+1)))
	# v = X[N - 1] - X[N]
	# v[N - 1] = 0.0
	# print(scipy.linalg.solve_triangular(np.transpose(C), v, lower=False))

	gradient = -np.array(gradient)
	ANS = -ANS
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
	N = 3

	INPUT = read_from_file("in.txt")

	X0 = np.array(INPUT[0])
	# X0 = initial_guess(N)
	# for i in range(N + 1):
	# 	for j in range(i):
	# 		# X0[i][j] = INPUT[0][i // 2][j // 2]
	# 		X0[i][j] = INPUT[0][i][j]
	# print(X0)

	# X0[1][0] = 1.0
	# X0[2][0] = 0.5
	# X0[2][1] = 1.0
	# X0[4][3] = 0.0

	modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
	for i in range(N + 1):
		for j in range(i):
			modifiable[i][j] = True
	# modifiable[1][0] = False
	# modifiable[2][0] = False
	# modifiable[2][1] = False
	# modifiable[4][3] = False

	Z = np.eye(N + 1, dtype='double')
	Z[0][0] = 0.0
	Z[1][0] = -1.0
	Z[1][1] = -1.0
	Z[2][0] = -1.0
	Z[2][1] = 1.0
	Z[2][2] = -1.0

	# S0 = np.ones((N), dtype='double')
	# S0 = []
	# modifiableS = []
	modifiableS = [True for i in range(N)]
	# modifiableS[0] = False
	INS = get_S(INPUT[1])
	S0 = [INS[i] for i in range(N)]
	
	fun = lambda v: calc(get_X(N, v, X0, S0))
	fun(get_v(N, X0, S0))
	# RES = scipy.optimize.minimize(fun, get_v(N, X0, S0), jac=True, options={"gtol": 0.0000000000000000000001}, method='BFGS')
	# RES = scipy.optimize.minimize(fun, get_v(N, X0, S0), jac=True, method='CG')
	# RES = scipy.optimize.minimize(fun, get_v(N, X0, S0), jac=True)

	# X1 = get_X(N, RES.x, X0, S0)[0]
	
	print("=================================================================")
	print(RES)

	# fun(RES.x)

	print(get_X(N, RES.x, X0, S0)[0])
	print(get_X(N, RES.x, X0, S0)[1])

	# with open("out.txt", "a") as f:
	# 	f.write(str(get_X(N, RES.x, X0, S0)[0]))
	# 	f.write("\n")
	# 	f.write(str(get_err2(get_X(N, RES.x, X0, S0)[1])))
	# 	f.write("\n\n\n")

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


