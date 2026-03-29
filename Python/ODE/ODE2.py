import numpy as np
import scipy
import math
import random

N = 4
V = np.zeros((N, N), dtype='double')
X = np.zeros((N, N), dtype='double')

S = 0.1
STEPS = 1000
G = []
W = []
M = np.zeros((N, N), dtype='double')

tab = []
e = np.zeros((N), dtype='double')
# def initialize_X():
	

def initialize():
	global N, V, X, G, W, M

	# X[1][0] = 1.0
	# X[2][0] = 0.3
	# X[2][1] = 0.8
	for i in range(N):
		for j in range(i):
			X[i][j] = random.random() - 1.0/2

	for i in range(N):
		for j in range(N):
			M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	
	V = np.linalg.cholesky(M)
	# print(V)
	# print(scipy.linalg.solve_triangular(V, np.identity(N), lower=True))

	M = np.dot(np.transpose(V), np.linalg.inv(M))

	G.append(np.zeros((N), dtype='double'))
	W.append(np.zeros((N), dtype='double'))

def next_step(k):
	global N, V, X, G, W, M, S, STEPS, tab, e
	h = 1.0*S/STEPS

	g = np.zeros((N), dtype='double')

	C = np.dot(X, W[k])
	if k % 100 == 0:
		tab.append(np.dot(e, W[k]))
	# if k % 100 == 0:
	# 	tab.append(-2*math.log(1 - (k*h)**2/2) - np.dot(W[k], W[k]))
	for i in range(N):
		g[i] = G[k][i] + h*math.exp(C[i] - np.dot(X[i], X[i])/2.0)
	G.append(g)
	# W.append(np.dot(M, z))
	W.append(scipy.linalg.solve_triangular(V, g, lower=True))

def steps():
	global STEPS, S, W
	for i in range(STEPS + 1):
		print(W[i])
		next_step(i)
	# print("RESULT:")
	# print(-2*math.log(1 - S*S/2))
	# print(-2*math.log(1 - S*S/2) - np.dot(W[STEPS], W[STEPS]))

initialize()
steps()

# for z in range(10000):
# 	print(z)

# 	tab.clear()
# 	for i in range(N):
# 		e[i] = random.random() - 1.0/2

# 	initialize()
# 	steps()

# 	for i in range(1, len(tab) - 1):
# 		if (tab[i] - tab[i - 1])*(tab[i + 1] - tab[i]) < -0.000001:
# 			print("BOOM!")
