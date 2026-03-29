# Covariance is exp(-||x - y||^2 / 2)
# Now we analyze the optimum of original target function

import numpy as np
import scipy
import random

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)

def K_function(x):
	D = len(x)
	# return [np.sin(j * x[j]) for j in range(D)]
	return [1.0 + np.sin(x[j] * (j+1)) * x[j] for j in range(D)]

def K(x, y):
	return np.exp(-np.dot(x - y, x - y) / 2.0)
	# return 4.0*np.exp(-np.dot(x - y, x - y) / 4.0) + np.exp(-np.dot(x - y, x - y))
	# return (1.0 + np.dot(x - y, x - y) / (2.0 * 1.0 * 1.0**2))**(-1)
	# return np.dot(K_function(x), K_function(y))

def get_b(X_initial, err2, X_evaluation, precision, X_target):
	X_initial = np.array(X_initial)
	N = len(X_initial)
	D = len(X_target)

	# indices 0 to N-1 are meant for initial data
	# index N is for evaluation point
	# index N+1 is for target point

	X = np.zeros((N + 2, D), dtype='double')
	for j in range(D):
		for i in range(N):
			X[i][j] = X_initial[i][j]
		X[N][j] = X_evaluation[j]
		X[N+1][j] = X_target[j]

	Sigma = np.zeros((N + 2, N + 2), dtype='double')
	for i in range(N + 2):
		for j in range(N + 2):
			Sigma[i][j] = K(X[i], X[j])
	for i in range(N):
		Sigma[i][i] += err2[i]
	Sigma[N][N] += 1.0/precision

	return scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+1][:(N+1)], assume_a='pos')


def G(b, X_initial, X_evaluation, X_intermediate, X_target):
	X_initial = np.array(X_initial)
	N = len(X_initial)

	Sigma_intermediate = np.zeros((N+1), dtype='double')
	for i in range(N):
		Sigma_intermediate[i] = K(X_initial[i], X_intermediate)
	Sigma_intermediate[N] = K(X_evaluation, X_intermediate)

	ANS = K(X_target, X_intermediate)
	ANS -= np.dot(Sigma_intermediate, b)

	return ANS

def target(X_initial, err2, X_evaluation, precision, X_target):
	b = get_b(X_initial, err2, X_evaluation, precision, X_target)

	W = np.zeros((D), dtype='double')
	for j in range(N):
		W += -(X_evaluation - X_initial[j]) * np.exp(-np.dot(X_evaluation - X_initial[j], X_evaluation - X_initial[j]) / 2.0) * b[j]
	W -= -(X_evaluation - X_target) * np.exp(-np.dot(X_evaluation - X_target, X_evaluation - X_target) / 2.0)
	LATENCY = np.dot(W, W) * np.exp(-np.dot(X_evaluation - X_target, X_evaluation - X_target) / 2.0)

	# print("LATENCY:", LATENCY)
	# LATENCY = 0.0

	return (G(b, X_initial, X_evaluation, X_target, X_target), LATENCY)


def printing(X_initial, err2, X_evaluation, precision, X_target, X_additional):
	b = get_b(X_initial, err2, X_evaluation, precision, X_target)
	original = G(b, X_initial, X_evaluation, X_evaluation, X_target)

	val = G(b, X_initial, X_evaluation, X_target, X_target)
	# print(original - val)
	if original - val < 0.0:
		return False
	
	# print("b:", b)
	# print("original:", original)
	for i in range(N):
		val = err2[i] * b[i]
		# print(original - val)
		# if original - val < 0.0:
		if original - val < 0.0:
			return False

	for k in range(len(X_additional)):
		val = G(b, X_initial, X_evaluation, X_additional[k], X_target)
		# print(original - val)
		# if original - abs(val) < 0.0:
		if original - val < 0.0:
			return False

	return True


def rand():
	return 1.0 - 2.0 * random.random()

def rand_err():
	# if random.random() < 0.5:
	# 	return np.exp(3.0 - 6.0*random.random())
	# return 0.0001
	return np.exp(3.0 - 20.0*random.random())

N = 2
M = 30
D = 1

for i in range(10000):
	print(i)

	X_initial = np.array([[rand() for j in range(D)] for i in range(N)])
	# X_initial = [[-0.85404753   ], [-0.08320190   ]]

	err2 = np.zeros((N), dtype='double')
	for i in range(N):
		err2[i] = rand_err()
		# err2[i] = 1.0
	# err2 = [5.89228798,      0.00010000     ]

	X_target = np.array([0.0 for j in range(D)])
	# X_target = [-0.79960979   ]

	# precision = 1.0/rand_err()
	# precision = 0.05390768904388641
	precision = 1.0

	x_additional = np.array([[rand() for j in range(D)] for i in range(M)])


	SUCCESS = False
	for k in range(100):
		v0 = np.array([rand() for j in range(D)])
		fun = lambda v: target(X_initial, err2, np.array(v), precision, X_target)[0]

		RES = scipy.optimize.minimize(fun, v0, options={"gtol": 0.0000000000000000000001})

		if printing(X_initial, err2, np.array(RES.x), precision, X_target, x_additional):
			SUCCESS = True
			break

	if not SUCCESS:
		print(RES)
		print("initial data:")
		print(X_initial)
		print("err2:")
		print(err2)
		print("evaluation point:")
		print(np.array(RES.x))
		print("precision:")
		print(precision)
		print("target point:")
		print(X_target)
		print("===================================")
		print(np.array(RES.x))
		(ANS, LATENCY) = target(X_initial, err2, np.array(RES.x), precision, X_target)
		print("LATENCY:", LATENCY)
		print("ANS:", ANS)
		printing(X_initial, err2, np.array(RES.x), precision, X_target, x_additional)

