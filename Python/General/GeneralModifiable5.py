# In this GeneralModifiable we assume that C[N][k] / X[N][k] = h = const
# However, the optimality conditions of a solution are different.

import numpy as np
import scipy
import random
import scipy.linalg
from pathlib import Path
 
np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f} \t'.format}, linewidth=200000, threshold=np.inf)

n = None
N = None
modifiableX = None
R2 = None


def get_v(X):
	global N, modifiableX

	v = []
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiableX[i][j]:
				v.append(X[i][j])

	return v


def get_X(v, X0):
	global N, modifiableX

	X = np.zeros((N + 1, N + 1), dtype='double')
	k = 0
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiableX[i][j]:
				X[i][j] = v[k]
				k += 1
			else:
				X[i][j] = X0[i][j]

	return X


# =============================================================================================================================

# Calculates value h and triple (X, sigma, Sigma, C)
def initialize(X, sigma):
	global n, N, P, R2, h

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
			Sigma[j][i] = Sigma[i][j]
	for i in range(N):
		Sigma[i][i] += sigma[i]
	C = scipy.linalg.cholesky(Sigma, lower=True)

	h = np.sqrt(R2 / (1.0 - np.dot(C[N, :n], C[N, :n])))
	print("factor h, must be at most 1.0:", h)

	return (X, sigma, Sigma, C)


# This calculates the error.
def calc(X, sigma):
	global n, N, P, R2, h

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
			Sigma[j][i] = Sigma[i][j]
	for i in range(N):
		Sigma[i][i] += sigma[i]
	C = scipy.linalg.cholesky(Sigma, lower=True)
	Cinv = np.linalg.inv(C)
	# print(C)
	
	p = np.zeros((N + 1), dtype='double')
	p[n:N] = 1.0 / sigma[n:]

	K = np.zeros((N + 1, N + 1), dtype='double')
	Kinv = np.zeros((N + 1, N + 1), dtype='double')
	K[n:, n:] = Sigma[n:, n:] - Sigma[n:, :n] @ np.linalg.inv(Sigma[:n, :n]) @ Sigma[:n, n:]
	Kinv[n:, n:] = np.linalg.inv(K[n:, n:])
	# print(K)
	# print(Kinv)

	# print(Kinv[N][N])
	print(Cinv[N, n:] / Cinv[N][N], Kinv[N, n:] / Kinv[N][N])
	# print(Kinv[N])

	m = np.sum(Cinv[N, n:N]) / (-np.sum(p[n:N]))
	p[N] = -1.0 / m
	# print(m)
	
	ERROR = 0.0
	for k in range(n, N):
		# print(m * p[k] + Cinv[N][k])
		ERROR += (m * p[k] + Cinv[N][k])**2
	
	for k in range(n, N):
		v = Sigma[(k + 1):N, (k + 1):] @ np.multiply(p[(k + 1):], X[(k + 1):, k])
		# print(v)
		ERROR += np.dot(v, v)

	ERROR += (np.dot(X[N][n:], X[N][n:]) - R2)**2

	print(ERROR)

	return ERROR



def iterate(X0, sigma0, Sigma0, C0):
	global n, N, P, R2, h

	# K = np.zeros((N + 1, N + 1), dtype='double')
	# for i in range(N + 1):
	# 	for j in range(N + 1):
	# 		K[i][j] = C0[N][max(i, j)] * X0[N][min(i, j)]
	# M = scipy.linalg.solve_triangular(C0.T, K.T, lower=False)
	# M = scipy.linalg.solve_triangular(C0.T, M.T, lower=False)
	
	# A = np.multiply(M, Sigma0)
	# U = np.flip(scipy.linalg.cholesky(-np.flip(A[n:N, n:N], (0, 1)), lower=True), (0, 1))

	e = np.zeros((N + 1), dtype='double')
	e[N] = 1.0
	a = scipy.linalg.solve_triangular(C0.T, e, lower=False)
	a *= C0[N][N]

	# K = Sigma0[n:(N + 1), n:(N + 1)] - Sigma0[n:(N + 1), 0:n] @ np.linalg.inv(Sigma0[0:n, 0:n]) @ Sigma0[0:n, n:(N + 1)]
	# print(a[:-1] - np.linalg.inv(K)[-1, :-1] / np.linalg.inv(K)[-1, -1])

	# print()
	# print("A = ")
	# print(A)
	# print("a = ")
	# print(a)
	# print(sigma0)




	X1 = np.copy(X0)
	sigma1 = np.copy(sigma0)
	Sigma1 = np.copy(Sigma0)
	C1 = np.copy(C0)

	sum_a = 0.0
	for k in range(n, N):
		sum_a += np.abs(a[k])
	for k in range(n, N):
		sigma1[k] = sum_a / (P * np.abs(a[k]))

	print("m =", sum_a/P)
	
	# for j in range(n):
	# 	# X1[n:N, j] = -np.linalg.solve(A[n:N, n:N], A[n:N, :n] @ X1[:n, j] + A[n:N, N] * X1[N][j])
	# 	X1[n:N, j] = scipy.linalg.solve_triangular(U, A[n:N, :n] @ X1[:n, j] + A[n:N, N] * X1[N][j], lower=False)
	# 	X1[n:N, j] = scipy.linalg.solve_triangular(U.T, X1[n:N, j], lower=True)

	for j in range(n, N - 1):
		# X1[N][j] = h * C0[N][j]
		# X1[(k + 1):N, k] = -np.linalg.solve(A[(k + 1):N, (k + 1):N], A[(k + 1):N, N] * X1[N][k])
		# X1[(j + 1):N, j] = scipy.linalg.solve_triangular(U[(j + 1 - n):, (j + 1 - n):], A[(j + 1):N, N] * X1[N][j], lower=False)
		# X1[(j + 1):N, j] = scipy.linalg.solve_triangular(U[(j + 1 - n):, (j + 1 - n):].T, X1[(j + 1):N, j], lower=True)
		X1[(j + 1):N, j] = np.divide(X1[N][j] * np.linalg.solve(Sigma0[(j + 1):N, (j + 1):N], Sigma0[(j + 1):N, N]), -a[(j + 1):N])
		# print(X1[(j + 1):N, j], np.divide(X1[N][j] * np.linalg.solve(Sigma0[(j + 1):N, (j + 1):N], Sigma0[(j + 1):N, N]), -a[(j + 1):N]))
		# print(np.multiply(-a[(j + 1):N], X1[(j + 1):N, j]) - X[N][j] * np.linalg.solve(Sigma0[(j + 1):N, (j + 1):N], Sigma0[(j + 1):N, N]))
	# X1[N][N - 1] = h * C0[N][N - 1]

	Sigma1 = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i + 1):
			Sigma1[i][j] = np.exp(-np.dot(X1[i] - X1[j], X1[i] - X1[j]) / 2.0)
			Sigma1[j][i] = Sigma1[i][j]
	for i in range(N):
		Sigma1[i][i] += sigma1[i]
	C1 = scipy.linalg.cholesky(Sigma1, lower=True)



	ANS = np.dot(X1[N, n:N], C1[N, n:N])
	# print("ANS:", ANS)

	return (X1, sigma1, Sigma1, C1)


# =============================================================================================================================

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


# =============================================================================================================================

X0 = None

# Reading from input file
# Indices 0 to M - 1 are meant for initial data.
# Indices M to N - 1 are meant for evaluation points.
# Index N is meant for output.
if False:
	INPUT = read_from_file("in.txt")
	X0 = np.array(INPUT[0])
	sigma2 = np.array(INPUT[1])
	N = len(X0) - 1

	n = None
	R2 = None

else:
	N = 4
	n = 0
	R2 = 0.1
	P = 1.0

	X0 = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i):
			X0[i][j] = 1.0/np.sqrt(j + 1)
			# X0[i][j] = random.random()
			# X0[i][j] = 0.0
			# X0[i][j] = 1.0/np.sqrt((j + 1) * (i + 1))
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
	# X0[N][N - 2] = random.random()
	# X0[N][N - 1] = random.random()
	# X0[N][N - 2] = 0.6
	# X0[N][N - 1] = 0.7
	# X0[N][0] = 1.04187201
	# X0[N][1] = 0.96769967
	# X0[N][2] = 0.48497041
	for j in range(N):
		X0[N][j] = random.random()
	# X0[N][N - 1] = np.sqrt(R2)

	sigma = np.zeros((N), dtype='double')
	# sigma[0] = 4.418927941
	# sigma[1] = 0.834729724
	# sigma[2] = 2.943724987
	sigma[0] = 1.0
	sigma[1] = 2.0
	sigma[2] = 3.0
	sigma[3] = 4.0
	# for i in range(n, N):
		# sigma[i] = (N - n) / P
		# sigma[i] = (N - n) / P * (1.0 + random.random())
	# for i in range(n, N):
	# 	sigma[i] = 10000000000000000.0
	# sigma[N - 1] = 1.0/P


# X = np.copy(X0)
# (X, sigma, Sigma, C) = initialize(X, sigma)

# for ITER in range(200):
# 	(X, sigma, Sigma, C) = iterate(X, sigma, Sigma, C)
	
# 	new_ANS = np.dot(X[N, n:N], C[N, n:N])


modifiableX = [[False for j in range(N + 1)] for i in range(N + 1)]
for i in range(n, N + 1):
	for j in range(n, i):
		modifiableX[i][j] = True

fun = lambda v: calc(get_X(v, X0), sigma)
RES = scipy.optimize.minimize(fun, get_v(X0), options={"gtol": 0.0000000000000000000001})

X = get_X(RES.x, X0)

print("DONE")
print(X)

with open("out.txt", "a") as f:
	f.write(str(X))
	f.write("\n")
	f.write(str(sigma))
	f.write("\n\n\n")
