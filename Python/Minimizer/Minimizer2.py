# A modified practical (still naive) implementation of Minimizer resulting from TrueGradient


import numpy as np
import random
from scipy.stats import ortho_group
from scipy.linalg import solve
import math

h = 1
factor = 1.618033988749894
# factor = 2.618033988749894

N = 2
A = np.zeros((N, N), dtype='double')

alpha = 1.0
times = []
grads = []
current_t = 0.0
current_x = np.zeros((N))
current_f = 0.0
comb = np.zeros((1))
start_index = 0
end_index = 1

diff_size = 0.0

def f(x):
	global A
	# v = np.dot(A, x)
	# return np.dot(v, v)
	# return math.sin(x[0])
	return (1 - x[0])**2 + 100*(x[1] - x[0]**2)**2

def df(x):
	global A
	# return 2*np.dot(np.transpose(A), np.dot(A, x))
	# return np.array([math.cos(x[0])])
	return np.array([-2*(1 - x[0]) - 400*x[0]*(x[1] - x[0]**2), 200*(x[1] - x[0]**2)])

def create_random_matrix():
	global A, N

	np.random.seed(seed=1)
	random.seed(1)
	A = ortho_group.rvs(dim = N)
	for i in range(N):
		c = random.random()
		c = c/(1 - c)
		for j in range(N):
			A[i][j] *= c

def sech2(d):
	if abs(d) > 15.0:
		return 0.0
	else:
		return 1.0/math.cosh(d)**2

def get_vector_comb():
	global comb, h, alpha, times, current_t, start_index, end_index

	end_index = len(times)
	start_index = end_index
	while start_index > 0 and end_index - start_index < 100 and (current_t - times[start_index - 1]) + h < 20.0:
		start_index -= 1

	K = end_index - start_index
	M = np.zeros((K, K))
	v = np.zeros((K))
	for i in range(start_index, end_index):
		for j in range(start_index, end_index):
			M[i - start_index][j - start_index] = sech2((times[i] - times[j])/alpha)
		v[i - start_index] = sech2((current_t - times[i])/alpha + h)
	comb = solve(M, v, assume_a='pos')


def iteration():
	global comb, N, h, alpha, factor, times, grads, current_f, current_t, current_x, start_index, end_index, diff_size
	
	alpha *= factor

	while True:
		get_vector_comb()
		# print(comb)

		diff = np.zeros(N)
		for i in range(start_index, end_index):
			diff += comb[i - start_index] * grads[i]
		diff *= h * alpha
		new_x = current_x - diff
		# print(math.sqrt(np.dot(current_x, current_x)), math.sqrt(np.dot(new_x, new_x)))
		new_f = f(new_x)
		diff_size = math.sqrt(np.dot(diff, diff))
		# print(current_x, new_x, diff_size)

		if new_f < current_f:
			current_t += h * alpha
			times.append(current_t)
			current_x = new_x
			current_f = new_f

			grads.append(df(current_x))
			break
		else:
			alpha /= factor**2
			continue


# create_random_matrix()
# print(A)
# A = np.identity(N)
# for i in range(N):
# 	A[i][i] = i + 1
current_x = np.zeros(N)
# current_x[0] = 1
current_t = 0.0
current_f = f(current_x)
times.append(0.0)
grads.append(df(current_x))

for i in range(1000):
	print(current_f, '\t', alpha, '\t', diff_size, '\t', end_index - start_index)
	iteration()


# x0 = np.zeros(N)
# x0[0] = 1
# x1 = np.zeros(N)
# x1[1] = 1
# print(f(x0))
# print(f(x1))
# print(f(x0 + x1))
