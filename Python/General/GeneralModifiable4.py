# In this GeneralModifiable we assume that C[N][k] / X[N][k] = h = const
# In fact: h = sqrt(1 - e^(-R2) rho2)/sqrt(R2)
# Here rho2 = 1.0 - C[N][0]**2 - ... - C[N][n-1]**2
# The algorithm is iterative!

import numpy as np
import scipy
import random
import scipy.linalg
from pathlib import Path

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.20f} \t'.format}, linewidth=200000, threshold=np.inf)

N = None
n = None
sigma2 = None
R2 = None


# =============================================================================================================================

# Initialize matrices Sigma, C
def initialize(X):
	global N, n, sigma2, R2

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
			Sigma[j][i] = Sigma[i][j]
		Sigma[i][i] += sigma2[i]
	C = scipy.linalg.cholesky(Sigma, lower=True)

	for j in range(n):
		Sigma[N][j] = np.exp(-(np.dot(X[N, :n] - X[j, :n], X[N, :n] - X[j, :n]) + R2) / 2.0)
		Sigma[j][N] = Sigma[N][j]
		C[N][j] = (Sigma[N][j] - np.dot(C[N, :j], C[j, :j])) / C[j][j]

	print("must be at most 1.0:", R2 + np.exp(-R2) * np.dot(C[N, :n], C[N, :n]))

	h = np.sqrt((1.0 - np.dot(C[N, :n], C[N, :n])) / R2)
	R2res = R2
	for j in range(n, N):
		Sigma[N][j] = np.exp(-(np.dot(X[N, :j] - X[j, :j], X[N, :j] - X[j, :j]) + R2res) / 2.0)
		Sigma[j][N] = Sigma[N][j]

		C[N][j] = (Sigma[N][j] - np.dot(C[N, :j], C[j, :j])) / C[j][j]
		X[N][j] = C[N][j] / h
		R2res -= X[N][j]**2
	
	W = np.zeros((N - n, N - n), dtype='double')
	for i in range(n, N):
		for j in range(n, N):
			W[i - n][j - n] = h * C[i][j] - Sigma[N][i] * X[i][j]

	w = np.zeros((N + 1), dtype='double')
	w[n:N] = np.linalg.solve(W.T, X[N, n:N])
	u = C.T @ w

	# Elements in index [N][N] should never be relevant in any way
	if R2res > 0.0:
		X[N][N] = np.sqrt(R2res)

	return (Sigma, C, u)


def iterate(X0, Sigma0, C0, u0):
	global N, n, sigma2, R2
	# We assume that X0 has been initialized i.e.
	# already has property C0[N][j] = h * X0[N][j].

	h = np.sqrt((1.0 - np.dot(C0[N, :n], C0[N, :n])) / R2)
	K = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			K[i][j] = C0[N][max(i, j)] * u0[min(i, j)]
	M = scipy.linalg.solve_triangular(C0.T, K.T, lower=False)
	M = scipy.linalg.solve_triangular(C0.T, M.T, lower=False)

	EPS = 0.000000000000000001
	for i in range(n, N):
		if M[i][i] > -EPS:
			M[i][i] = -EPS
			print("bum")
	
	A = np.multiply(M, Sigma0)
	U = np.flip(scipy.linalg.cholesky(-np.flip(A[n:N, n:N], (0, 1)), lower=True), (0, 1))

	X1 = np.copy(X0)
	Sigma1 = np.copy(Sigma0)
	C1 = np.copy(C0)
	
	for j in range(n):
		# X1[n:N, j] = -np.linalg.solve(A[n:N, n:N], A[n:N, :n] @ X1[:n, j] + A[n:N, N] * X1[N][j])
		X1[n:N, j] = scipy.linalg.solve_triangular(U, A[n:N, :n] @ X1[:n, j] + A[n:N, N] * X1[N][j], lower=False)
		X1[n:N, j] = scipy.linalg.solve_triangular(U.T, X1[n:N, j], lower=True)


	# Now the main loop
	R2res = R2
	for k in range(n, N):

		# First update Sigma1[k][:k] and C1[k][:k]
		for j in range(k):
			Sigma1[k][j] = np.exp(-np.dot(X1[k] - X1[j], X1[k] - X1[j]) / 2.0)
			Sigma1[j][k] = Sigma1[k][j]
			C1[k][j] = (Sigma1[k][j] - np.dot(C1[k, :j], C1[j, :j])) / C1[j][j]

		# Next update C1[k][k] and Sigma1[N][k] and C1[N][k]
		C1[k][k] = np.sqrt(max(Sigma1[k][k] - np.dot(C1[k, :k], C1[k, :k]), sigma2[k]))
		Sigma1[N][k] = np.exp(-(np.dot(X1[N, :k] - X1[k, :k], X1[N, :k] - X1[k, :k]) + R2res) / 2.0)
		Sigma1[k][N] = Sigma1[N][k]
		C1[N][k] = (Sigma1[N][k] - np.dot(C1[N, :k], C1[k, :k])) / C1[k][k]

		# Now that we have C1[N][k], we can update X1[N][k] and R2res
		X1[N][k] = C1[N][k] / h
		R2res -= X1[N][k]**2

		# Finally, we update X1[(k + 1):N, k].
		if k < N - 1:
			# X1[(k + 1):N, k] = -np.linalg.solve(A[(k + 1):N, (k + 1):N], A[(k + 1):N, N] * X1[N][k])
			X1[(k + 1):N, k] = scipy.linalg.solve_triangular(U[(k + 1 - n):, (k + 1 - n):], A[(k + 1):N, N] * X1[N][k], lower=False)
			X1[(k + 1):N, k] = scipy.linalg.solve_triangular(U[(k + 1 - n):, (k + 1 - n):].T, X1[(k + 1):N, k], lower=True)

	u1 = ((A[N] @ X1) + X1[N]) / h
	u1[N] = 0.0


	if R2res > 0.0:
		X1[N][N] = np.sqrt(R2res)
	ANS = R2res / 2.0
	print("ANS:", ANS)

	return (X1, Sigma1, C1, u1)


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
	N = 200
	n = 3
	R2 = 0.01
	P = 100.0

	X0 = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i):
			X0[i][j] = 1.0/np.sqrt(j + 1)
			# X0[i][j] = 1.0/np.sqrt((j + 1) * (i + 1))
	X0[0][0] = 0.5798527
	X0[0][1] = 0.32425464266
	X0[0][2] = 0.163347653754
	X0[1][0] = -0.57285235
	X0[1][1] = -0.067735435654
	X0[1][2] = -0.35463773
	X0[2][0] = 0.867143785681
	X0[2][1] = 0.538572989
	X0[2][2] = -0.14325635
	X0[N][0] = 0.0
	X0[N][1] = 0.0
	X0[N][2] = 0.0
	# X0[N][0] = 1.04187201
	# X0[N][1] = 0.96769967
	# X0[N][2] = 0.48497041

	sigma2 = np.zeros((N + 1), dtype='double')
	sigma2[0] = 4.418927941
	sigma2[1] = 0.834729724
	sigma2[2] = 2.943724987
	for i in range(n, N):
		sigma2[i] = (N - n) / P


X = np.copy(X0)
(Sigma, C, u) = initialize(X)

EPS = 0.0000000001
BEST = 100000.0
SECOND_BEST = 1000000.0
ITER = 0

while ITER < 20:
	(X, Sigma, C, u) = iterate(X, Sigma, C, u)
	
	new_ANS = X[N][N]**2

	if new_ANS < BEST:
		SECOND_BEST = BEST
		BEST = new_ANS
	elif new_ANS < SECOND_BEST:
		SECOND_BEST = new_ANS
	else:
		ITER += 1
	# if  SECOND_BEST - BEST < EPS:
	# 	ITER += 1
	# else:
	# 	ITER = 0

print("DONE")

with open("out.txt", "a") as f:
	f.write(str(X))
	f.write("\n")
	f.write(str(sigma2))
	f.write("\n\n\n")
