# Some funky business with kernels. Nothing interesting tbh

import numpy as np
import random
from scipy.linalg import solve, solve_triangular
import math


def kernel(d):
	return 1.0/math.cosh(d)**2
	# return math.exp(-d*d)

K = 4
MAX = 2
h = 0.9

times = []
for i in range(K):
	times.append(random.random()*MAX)
times.sort()
times.reverse()

# times = [MAX - 2*h, MAX - 3*h, MAX - 5*h, MAX - 7*h]
# times = [MAX - h * (1 + k) for k in range(K)]

print(times)


def vec(t):
	x = math.tanh(t)
	return np.array([1, math.sqrt(2) * x, math.sqrt(3) * x**2, 2 * x**3, math.sqrt(5) * x**4, math.sqrt(6) * x**5]) / math.cosh(t)**2

def kernel2(t1, t2):
	return np.dot(vec(t1), vec(t2))


M = np.zeros((K, K))
v = np.zeros((K))
for i in range(K):
	for j in range(K):
		M[i][j] = kernel(times[i] - times[j])
		# M[i][j] = kernel2(times[i], times[j])
	v[i] = kernel(MAX - times[i])
	# v[i] = kernel2(MAX, times[i])
comb = solve(M, v, assume_a='pos')
# comb = solve(M, v)

print(M)
print(np.linalg.inv(M))


def f(t):
	return math.tanh(t)
	# return t

U = np.zeros((K, K))
for k in range(K):
	tab = [[0.0 for j in range(K)] for i in range(K)]

	for i in range(K):
		if i == k:
			tab[i][i] = 1.0
		else:
			tab[i][i] = 0.0
	
	for d in range(1, K):
		for i in range(K - d):
			j = i + d
			tab[i][j] = (tab[i + 1][j] - tab[i][j - 1])/(f(times[j]) - f(times[i]))

	for i in range(K):
		U[k][i] = tab[0][i]
			
# print(U)

UMU = np.dot(np.dot(np.transpose(U), M), U)
print(UMU)
comb2 = solve(UMU, np.dot(np.transpose(U), v), assume_a='pos')
print(np.linalg.inv(UMU))

# print(comb)
# print(solve_triangular(U, comb, lower=False))	

# print(np.dot(vec(0), (vec(0.02) - 2*vec(0.01) + vec(0)))/2.0)
# print(vec(0.03))
# print(vec(0.02) * comb2[0] + 100*(vec(0.02) - vec(0.01)) * comb2[1] + 5000 * (vec(0.02) - 2*vec(0.01) + vec(0)) * comb2[2])
# print(vec(0.02) * comb2[0] + 100*(vec(0.02) - vec(0.01)) * 0.01 + 5000 * (vec(0.02) - 2*vec(0.01) + vec(0)) * 0.01**2 / 2)
print(comb2)
FACTOR = 1.0
for i in range(K):
	print(comb2[i] * FACTOR)
	FACTOR /= (MAX - times[i])
