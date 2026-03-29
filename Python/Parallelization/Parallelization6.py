# Covariance is exp(-||x - y||^2 / 2)
# Now we analyze the extrema of G function (named calc here)

import numpy as np
import scipy
import random

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)


# X_initial, err2, precision, X_target are fixed
# X_evaluation, X_intermediate are variables
def calc(X_initial, err2, X_evaluation, precision, X_intermediate, X_target):
	X_initial = np.array(X_initial)
	N = len(X_initial)
	D = len(X_target)

	# indices 0 to N-1 are meant for initial data
	# index N is for evaluation point
	# index N+1 is for intermediate point
	# index N+2 is for target point

	X = np.zeros((N + 3, D), dtype='double')
	for j in range(D):
		for i in range(N):
			X[i][j] = X_initial[i][j]
		X[N][j] = X_evaluation[j]
		X[N+1][j] = X_intermediate[j]
		X[N+2][j] = X_target[j]

	Sigma = np.zeros((N + 3, N + 3), dtype='double')
	for i in range(N + 3):
		for j in range(N + 3):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	for i in range(N):
		Sigma[i][i] += err2[i]
	Sigma[N][N] += 1.0/precision	

	ANS = Sigma[N+1][N+2] - np.dot(Sigma[N+1][:(N+1)], scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+2][:(N+1)], assume_a='pos'))

	a = scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+1][:(N+1)], assume_a='pos')
	b = scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+2][:(N+1)], assume_a='pos')


	LATENCY = 0.0
	W = np.zeros((D), dtype='double')
	for j in range(N + 1):
		W += -(X[N+1] - X[j]) * Sigma[N+1][j] * b[j]
	W -= -(X[N+1] - X[N+2]) * Sigma[N+1][N+2]
	LATENCY -= np.dot(W, W) * np.exp(1.2 * np.dot(X[N+1] - X[N+2], X[N+1] - X[N+2]))

	# W = np.zeros((D), dtype='double')
	# for j in range(N + 1):
	# 	W += -(X[N] - X[j]) * Sigma[N+1][j] * b[j]
	# W -= -(X[N+1] - X[N+2]) * Sigma[N+1][N+2]
	# LATENCY -= np.dot(W, W) * np.exp(1.2 * np.dot(X[N+1] - X[N+2], X[N+1] - X[N+2]))

	return (ANS, LATENCY)


def rand():
	return 2.0 - random.random()*4.0

def rand_err():
	return np.exp(3.0 - 6.0*random.random())


N = 0
D = 2

for i in range(1):

	X_initial = np.array([[rand() for j in range(D)] for i in range(N)])
	err2 = np.zeros((N), dtype='double')
	for i in range(N):
		err2[i] = rand_err()
	X_target = np.array([0.0 for j in range(D)])
	precision = 1.0/rand_err()

	v0 = np.array([rand() for j in range(2*D)])

	fun = lambda v: -calc(X_initial, err2, np.array(v[:D]), precision, np.array(v[D:]), X_target)[0]

	RES = scipy.optimize.minimize(fun, v0, options={"gtol": 0.0000000000000000000001})
	print(RES)
	print("initial data:")
	print(X_initial)
	print("err2:")
	print(err2)
	print("evaluation point:")
	print(np.array(RES.x[:D]))
	print("precision:")
	print(precision)
	print("intermediate point:")
	print(np.array(RES.x[D:]))
	print("target point:")
	print(X_target)
	print("===================================")
	print(np.array([RES.x[:D], RES.x[D:]]))
	(ANS, LATENCY) = calc(X_initial, err2, np.array(RES.x[:D]), precision, np.array(RES.x[D:]), X_target)
	print(ANS)
	print(LATENCY)

