# Covariance is exp(-||x - y||^2 / 2)
# Now we analyze the optimum of original target function
# Checking the conditions on finite sets (arbitrary covariance matrices!)

import numpy as np
import scipy
import random
from pprint import pprint
import scipy.linalg



np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

def K_function(x):
	D = len(x)
	# return [np.sin(j * x[j]) for j in range(D)]
	return [1.0 + np.sin(x[j] * (j+1)) * x[j] for j in range(D)]

def K_gen(x, y):
	return np.exp(-np.dot(x - y, x - y) / 2.0)
	# return 4.0*np.exp(-np.dot(x - y, x - y) / 4.0) + np.exp(-np.dot(x - y, x - y))
	# return (1.0 + np.dot(x - y, x - y) / (2.0 * 1.0 * 1.0**2))**(-1)
	# return (np.dot(x, y) + 1.0)**1
	# return np.dot(K_function(x), K_function(y))

def generate_covariance_matrix(n):
	global C, D

	C = np.array([[random.normalvariate(0.0, 1.0) for j in range(D)] for i in range(n)])
	# C = np.array([[abs(random.normalvariate(0.0, 1.0)) for j in range(D)] for i in range(n)])
	# C[n-1] = [0.0, 1.0]

	# for i in range(n):
	# 	l = np.sqrt(np.dot(C[i], C[i]))
	# 	C[i] /= l

	# 	max_element = C[i][i]
	# 	max_index = i
	# 	for j in range(n):
	# 		if C[i][j] > max_element:
	# 			max_element = C[i][j]
	# 			max_index = j
	# 	C[i][max_index] = C[i][i]
	# 	C[i][i] = max_element

	# F = np.array([[random.normalvariate(0.0, 1.0) for j in range(n)] for i in range(2*n)])
	# W = np.array([abs(random.normalvariate(0.0, 1.0)) for i in range(2*n)])

	# G = np.array([[0.0 for j in range(n)] for i in range(n)])
	# for i in range(n):
	# 	for j in range(n):
	# 		for k in range(2*n):
	# 			G[i][j] += W[k] * np.cos(np.dot(C[i] - C[j], F[k]))

	return np.dot(C, np.transpose(C))
	# return np.array([[K_gen(C[i], C[j]) for j in range(n)] for i in range(n)])
	# return G

def get_b(K, err2, precision, index_evaluation, N, M):

	# indices 0 to N-1 are initial data
	# index N is the evaluation point
	# index N+1 is the target point

	K2 = np.zeros((N + 2, N + 2), dtype='double')
	for i in range(N):
		for j in range(N):
			K2[i][j] = K[i][j]
		K2[i][i] += err2[i]
	for i in range(N):
		K2[N][i] = K[index_evaluation][i]
		K2[i][N] = K[i][index_evaluation]
		K2[N+1][i] = K[N+M][i]
		K2[i][N+1] = K[i][N+M]
	K2[N][N] = K[index_evaluation][index_evaluation] + 1.0/precision
	K2[N][N+1] = K[index_evaluation][N+M]
	K2[N+1][N] = K[N+M][index_evaluation]
	K2[N+1][N+1] = K[N+M][N+M]
	
	return scipy.linalg.solve(K2[:(N+1), :(N+1)], K2[N+1][:(N+1)], assume_a='pos')


def G(K, b, index_evaluation, index_intermediate, N, M):
	ANS = K[index_intermediate][N+M]
	for i in range(N):
		ANS -= K[i, index_intermediate] * b[i]
	ANS -= K[index_evaluation][index_intermediate] * b[N]

	return ANS


def calc(K, err2, precision, index_evaluation, N, M, printing):
	b = get_b(K, err2, precision, index_evaluation, N, M)
	value_at_evaluation = abs(G(K, b, index_evaluation, index_evaluation, N, M))
	# value_at_evaluation = G(K, b, index_evaluation, index_evaluation, N, M)
	target = G(K, b, index_evaluation, N+M, N, M)

	for i in range(N + M + 1):
		value_at_i = abs(G(K, b, index_evaluation, i, N, M))
		# value_at_i = G(K, b, index_evaluation, i, N, M)
		if printing:
			print("i =", i, "\tG =", value_at_i, "\toriginal =", value_at_evaluation)
		if value_at_evaluation < value_at_i:
			return (False, target)

	return (True, target)


def rand():
	return 1.0 - 2.0 * random.random()

def rand_err():
	return np.exp(3.0 - 6.0*random.random())


# indices 0 to N-1 are initial data
# indices N to N+M-1 are additional points
# index N+M is the target 
N = 0
M = 2
D = 2
C = None

for i in range(100000):
	print(i)

	K = generate_covariance_matrix(N + M + 1)

	err2 = np.zeros((N), dtype='double')
	for i in range(N):
		# err2[i] = rand_err()
		err2[i] = 1.0
	# precision = 1.0/rand_err()
	precision = 1.0



	best_index = N+M
	best_target = K[N+M][N+M]
	all_results = [calc(K, err2, precision, k, N, M, False) for k in range(N+M+1)]

	lowest_at = N+M
	lowest_value = K[N+M][N+M]
	working_amt = 0

	for k in range(N+M+1):
		if all_results[k][1] < lowest_value:
			lowest_value = all_results[k][1]
			lowest_at = k
		if all_results[k][0]:
			working_amt += 1

	# if working_amt >= 2:
	if (working_amt >= 1 and not all_results[lowest_at][0]) or (working_amt >= 2):
	# if working_amt != 1 or not all_results[lowest_at][0]:
		print("WRONG")
		print(C)
		pprint(all_results)
		print(lowest_value)
		print(err2)
		print(precision)
		print(K)
		print("========= PRINTING")
		for k in range(N+M+1):
			print("k =", k)
			calc(K, err2, precision, k, N, M, True)
		break
	
