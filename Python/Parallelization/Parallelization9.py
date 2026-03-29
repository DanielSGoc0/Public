# Covariance is exp(-||x - y||^2 / 2)
# Looking for the "right" symmetric matrices.


import numpy as np
import scipy
import random
from pprint import pprint
import scipy.linalg
import copy
import sys



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

	return np.dot(C, np.transpose(C))**3
	# return np.array([[K_gen(C[i], C[j]) for j in range(n)] for i in range(n)])
	# return G


def check(K, err2):

	WORKS = True
	for index_target in range(N):

		# Firstly we find the index_evaluation.
		# It is the maximizer of K[target][eval]^2 / (K[eval][eval] + err2)
		index_eval = index_target
		optimum = K[index_target][index_target]**2 / (K[index_target][index_target] + err2)

		for i in range(N):
			val = K[index_target][i]**2 / (K[i][i] + err2)
			if val > optimum:
				optimum = val
				index_eval = i
		
		# Now check the condition.
		bound = abs(K[index_eval][index_target]) * err2 / (K[index_eval][index_eval] + err2)

		for j in range(N):
			if j == index_eval:
				continue
		
			val  = abs(K[j][index_target] - K[j][index_eval] * K[index_eval][index_target] / (K[index_eval][index_eval] + err2))

			if val > bound:
				WORKS = False
				break

		if not WORKS:
			break

	return WORKS

def check2(K):

	WORKS = True
	for power in range(-100, 100):
		err2 = 1.1**power

		if not check(K, err2):
			WORKS = False
			break

	return WORKS

def dominant(K):

	for i in range(N):
		for j in range(N):
			if i == j:
				continue
			if abs(K[i][j]) > K[i][i]:
				return False
	return True



def all_eps(N):
	res = []
	current = [-1 for _ in range(N)]

	for _ in range(3**N - 1):
		res.append(copy.deepcopy(current))

		for k in range(N):
			if current[k] == -1:
				current[k] = 1
				break
			elif current[k] == 1:
				current[k] = 0
				break
			else:
				current[k] = -1

	return res

def find_optimal(K, target, precision):
	N = len(K)

	answers = []

	for eps in all_eps(N):

		S = 0
		for el in eps:
			S += el*el

		a = np.zeros((S), dtype='double')
		A = np.zeros((S, S), dtype='double')
		e = np.zeros((S), dtype='double')

		counteri = 0
		for i in range(N):
			if eps[i] == 0:
				continue

			a[counteri] = K[i][target]
			e[counteri] = eps[i]

			counterj = 0
			for j in range(N):
				if eps[j] == 0:
					continue

				A[counteri][counterj] = K[i][j]

				counterj += 1
			counteri += 1

		Ainv = np.linalg.solve(A, np.eye(S))
		M = np.dot(e, np.dot(Ainv, a)) / (precision + np.dot(e, np.dot(Ainv, e)))

		# A sanity check
		if M <= 0.0:
			continue

		WORKS = True

		v = np.dot(Ainv, a) - M * np.dot(Ainv, e)
		# Check that all p_k are non-negative
		for k in range(S):
			if e[k] * v[k] < 0.0:
				WORKS = False
				break
		if not WORKS:
			continue

		# And now the greatest test:
		for i in range(N):
			if eps[i] != 0:
				continue

			SUM = K[i][target]
			counterj = 0
			for j in range(N):
				if eps[j] == 0:
					continue
				SUM -= K[i][j] * v[counterj]

				counterj += 1

			if abs(SUM) > M:
				WORKS = False
				break

		if WORKS:
			answers.append([eps, v, M])

			for i in range(N):
				if eps[i] != 0:
					continue

				SUM = K[i][target]
				counterj = 0
				for j in range(N):
					if eps[j] == 0:
						continue
					SUM -= K[i][j] * v[counterj]

					counterj += 1

				print(SUM, M)
				if abs(SUM) > M:
					WORKS = False
					break

	return answers


def rand():
	return 1.0 - 2.0 * random.random()

def rand_err():
	return np.exp(3.0 - 6.0*random.random())


# indices 0 to N-1 are initial data
# indices N to N+M-1 are additional points
# index N+M is the target 
N = 3
D = 3
C = None

for i in range(1):
	print(i)

	K = generate_covariance_matrix(N)
	# precision = rand_err()

	# err2 = rand_err()
	# if check(K, err2):
	# 	print("===========")
	# 	print(K)
	# 	print(err2)

	# if check2(K) != dominant(K):
	# 	print("===========")
	# 	print(K)

	# if len(find_optimal(K, 0, precision)) != 1:
	# 	print("===========")
	# 	print(K)
	# 	print(precision)
	# 	break

	# K = [[0.05429126, -0.19052639, 0.19515997], [-0.19052639, 0.80027209, -0.95591485], [0.19515997, -0.95591485, 1.64107734]]
	# K = [[0.01, -0.01, -0.01], [-0.01, 1.0, 0.01], [-0.01, 1.0, 100.0]]
	# K = [[0.01362774, -0.02715800, -0.06765395], [-0.02715800, 0.08248598, 0.36064315], [-0.06765395, 0.36064315, 91.46188479]]
	# K = [[0.01, -0.03, -0.07], [-0.03, 0.1, 0.4], [-0.07, 0.4, 15.0]]
	# K = [[0.01, -0.03, -0.07], [-0.03, 0.1, 0.4], [-0.07, 0.4, 15.0]]
	# print(np.linalg.eigh(K))

	l = []
	for k in range(2):
	# for k in range(10, 11):
		precision = 10.0**(2*k - 1)
		# print("PRECISION =", precision)

		answers = find_optimal(K, 0, precision)

		if len(answers) != 1:
			print(K)
			sys.exit("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		
		for (eps, v, M) in answers:
			print(np.array(eps))
			print(np.array(v))
			print(M)
			# print(len(v))
			l.append(len(v))

	if l[0] > l[1]:
		print(np.array(K))
		break
