# Calculates the asymptotic limits of determinants, as the distance between points goes to zero.

import numpy as np
import random
import scipy
import itertools

import scipy.special

N = 4
D = 2

MIN = 100000
MAX = 0

def area(v0, v1, v2):
	return np.dot(v0 - v1, v0 - v1) * np.dot(v0 - v2, v0 - v2) - np.dot(v0 - v1, v0 - v2)**2

def custom_sum(X, s):
	global N

	S = 0.0
	for perm in itertools.permutations(s):
		Sigma = np.zeros((N, N), dtype = 'double')

		for i in range(N):
			for j in range(N):
				Sigma[i][j] = (-np.dot(X[i] - X[j], X[i] - X[j])/2)**perm[i] / scipy.special.factorial(perm[i])

		S += np.linalg.det(Sigma)
		
	return S

def custom_sum2(X, s):
	global N

	S = 0.0
	for perm in itertools.permutations(s):
		Sigma = np.zeros((N, N), dtype = 'double')

		for i in range(N):
			for j in range(N):
				Sigma[i][j] = (np.dot(X[j], X[i]))**perm[i] / scipy.special.factorial(perm[i])

		S += np.linalg.det(Sigma)
		
	return S


for k in range(1):
	X = np.array([[random.random() for j in range(D)] for i in range(N)])
	X[0] = [0, 1]
	X[1] = [1, 0]
	X[2] = [1, 1]
	X[3] = [0, 0]
	# X[N - 1] = 0
	B = np.zeros((N, N), dtype='double')

	# w0 = custom_sum(X, [0, 1, 1, 2]) / 2
	# print(w0)
	# w1 = custom_sum(X, [1, 1, 1, 1]) / 24
	# print(w1)
	# print(w0 + w1)
	# print(area(X[0], X[1], X[2]))
	print(np.dot(X[0] - X[3], X[0] - X[3]) * area(X[3], X[1], X[2]) / 2)
	print()

	# print(area(X[1], X[2], X[3]))
	

	for s in range(10):
		# A = np.zeros((N - 1, N - 1), dtype='double')

		Sigma = np.zeros((N, N), dtype = 'double')

		for i in range(N):
			for j in range(N):
				Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)

		# for i in range(N - 1):
		# 	for j in range(N - 1):
		# 		A[i][j] = np.dot(X[i] - X[N - 1], X[j] - X[N - 1])

		detSigma = np.linalg.det(Sigma)
		# detA = np.linalg.det(A)

		# # v = detSigma / detA
		v = detSigma
		# print(v)
		# MIN = min(MIN, v)
		# MAX = max(MAX, v)

		# w0 = custom_sum2(X, [0, 1, 1, 2]) / 2
		# # # print(w0)
		# w1 = custom_sum2(X, [1, 1, 1, 1]) / 24
		# # # print(w1)
		# v = w0 + w1

		print(v * 2**(14*s))

		# for j in range(N):
		X[0] /= 2.0
		X[1] /= 4.0
		X[2] /= 8.0


# print("MAX: ", MAX)
# print("MIN: ", MIN)
	
