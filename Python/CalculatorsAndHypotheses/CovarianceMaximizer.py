# Copied from GeneralModifiable4.py
# A method that takes the output point to evolve as in GeneralModifiable4,
# but the evaluation point is always a maximizer of expression:
# K(X_eval, X_out)^2 / (K(X_eval, X_eval) + 1/sigma_eval^2)
import numpy as np
import scipy
import random
import scipy.linalg
from pathlib import Path

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.12f} \t'.format}, linewidth=200000, threshold=np.inf)

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
			Sigma[i][j] = np.exp(-np.dot(X[i, :n] - X[j, :n], X[i, :n] - X[j, :n]) / 2.0)
			Sigma[j][i] = Sigma[i][j]
		Sigma[i][i] += sigma2[i]
	C = np.zeros((N + 1, N + 1), dtype='double')
	C[:n, :n] = scipy.linalg.cholesky(Sigma[:n, :n], lower=True)

	print("must be at most 1.0:", R2 + np.exp(-R2) * np.dot(C[N, :n], C[N, :n]))
	
	for j in range(n):
		Sigma[N][j] = np.exp(-(np.dot(X[N, :n] - X[j, :n], X[N, :n] - X[j, :n]) + R2) / 2.0)
		Sigma[j][N] = Sigma[N][j]
		C[N][j] = (Sigma[N][j] - np.dot(C[N, :j], C[j, :j])) / C[j][j]

	return (Sigma, C)


# Calculates the expression K(x_eval, x_out)^2 / (K(x_eval, x_eval) + 1 / sigma_eval^2)
def expression(X, Sigma, C, x_eval, k):
	global N, n, sigma2, R2

	x_eval = np.array(x_eval)
	Sigma_eval = np.zeros((k), dtype='double')
	for i in range(k):
		Sigma_eval[i] = np.exp(-np.dot(x_eval - X[i, :k], x_eval - X[i, :k]) / 2.0)
	Sigma_evalout = np.exp(-(np.dot(x_eval - X[N, :k], x_eval - X[N, :k]) + R2 - np.dot(X[N, n:k], X[N, n:k])) / 2.0)

	C_eval = scipy.linalg.solve_triangular(C[:k, :k], Sigma_eval, lower=True)

	return (Sigma_evalout - np.dot(C_eval, C[N, :k]))**2 / (1.0 + sigma2[k] - np.dot(C_eval, C_eval))
	# return Sigma_evalout - np.dot(C_eval, C[N, :k])
	# return (Sigma_evalout - np.dot(C_eval, C[N, :k])) / np.sqrt(1.0 + sigma2[k] - np.dot(C_eval, C_eval))
	# return np.dot(C_eval, C_eval)


def iterate(X, Sigma, C, k):
	global N, n, sigma2, R2
	# We assume that X0 satisfies C0[N][j] = h * X0[N][j] for all n <= j < k.
	# Our goal is to find the k-th evaluation point.
	# Moreover, we know all the entries X0[:k], Sigma[:k, :k], C0[:k]
	# as well as X0[N, :k], Sigma[N, :k], C0[N, :k]

	h = np.sqrt((1.0 - np.dot(C[N, :n], C[N, :n])) / R2)

	expr = lambda v: -expression(X, Sigma, C, v, k)
	RES = scipy.optimize.minimize(expr, X[N, :k], jac=False, options={"gtol": 0.0000000000000000000001}, method='BFGS')
	# RES = scipy.optimize.minimize(expr, np.ones((k), dtype='double'), jac=False, options={"gtol": 0.0000000000000000000001}, method='BFGS')

	# The result RES.x is the next evaluation point
	X[k, :k] = np.array(RES.x)

	# Now we have to update the relevant fields in Sigma and C
	for j in range(k, N + 1):
		X[k][j] = 0.0

	for i in range(k):
		Sigma[k][i] = np.exp(-np.dot(X[k] - X[i], X[k] - X[i]) / 2.0)
		Sigma[i][k] = Sigma[k][i]
	Sigma[k][k] = 1.0 + sigma2[k]
	Sigma[N][k] = np.exp(-(np.dot(X[N, :k] - X[k, :k], X[N, :k] - X[k, :k]) + R2 - np.dot(X[N, n:k], X[N, n:k])) / 2.0)
	Sigma[k][N] = Sigma[N][k]

	for j in range(k):
		C[k][j] = (Sigma[k][j] - np.dot(C[k, :j], C[j, :j])) / C[j][j]
	C[k][k] = np.sqrt(max(Sigma[k][k] - np.dot(C[k, :k], C[k, :k]), sigma2[k]))
	C[N][k] = (Sigma[N][k] - np.dot(C[N, :k], C[k, :k])) / C[k][k]

	# And finally, update X[N][k]
	X[N][k] = C[N][k] / h


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
	N = 100
	n = 3
	R2 = 0.5
	P = 10.0

	X0 = np.zeros((N + 1, N + 1), dtype='double')
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
		# sigma2[i] = 0.0

X = np.copy(X0)
(Sigma, C) = initialize(X)
for k in range(n, N):
	print(k)
	iterate(X, Sigma, C, k)
X[N][N] = np.sqrt(max(R2 - np.dot(X[N, n:N], X[N, n:N]), 0.0))

# print(X)

with open("out.txt", "a") as f:
	f.write(str(X))
	f.write("\n")
	f.write(str(sigma2))
	f.write("\n\n\n")
