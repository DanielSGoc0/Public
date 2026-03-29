import numpy as np
import scipy
import math

N = 3
V = np.zeros((N, N), dtype='double')
X = np.zeros((N, N), dtype='double')

T = 1
STEPS = 1000
Z = []
W = []
M = np.zeros((N, N), dtype='double')

def initialize():
	global N, V, X, Z, W, M

	X[1][0] = 0.4385
	X[2][0] = 0.6541
	X[2][1] = 0.2786

	for i in range(N):
		for j in range(N):
			M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	
	V = np.linalg.cholesky(M)
	# print(V)

	M = np.dot(np.transpose(V), np.linalg.inv(M))

	Z.append(np.zeros((N), dtype='double'))
	W.append(np.zeros((N), dtype='double'))

def next_step(k):
	global N, V, X, Z, W, M, T, STEPS
	h = 1.0*T/STEPS
	t = (2.0*k + 1)*T/(2*STEPS)
	factor = h * math.cosh(t/math.sqrt(2))**(-2)

	z = np.zeros((N), dtype='double')

	C = np.dot(X, W[k])
	for i in range(N):
		z[i] = Z[k][i] + math.exp(C[i] - np.dot(X[i], X[i])/2.0) * factor
	Z.append(z)
	# W.append(np.dot(M, z))
	W.append(scipy.linalg.solve_triangular(V, z, lower=True))


def steps():
	global STEPS, T, W
	for i in range(STEPS):
		# print(W[i])
		next_step(i)
	print("RESULT:")
	print(4*math.log(math.cosh(T/math.sqrt(2))))
	print(4*math.log(math.cosh(T/math.sqrt(2))) - np.dot(W[STEPS], W[STEPS]))

initialize()
steps()
