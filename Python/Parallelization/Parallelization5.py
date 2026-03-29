# Covariance is exp(-||x - y||^2 / 2)
# This one is for analyzing the G function.

import numpy as np
import scipy
import random

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=20000)


def target(X_initial, err2, X_evaluation, precision, X_target):
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
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	for i in range(N):
		Sigma[i][i] += err2[i]
	Sigma[N][N] += 1.0/precision

	return np.dot(Sigma[N+1][:(N+1)], scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+1][:(N+1)], assume_a='pos'))


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

	# conditioned_variance = Sigma[N+1][N+1] - np.dot(Sigma[N+1][:(N+1)], scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+1][:(N+1)], assume_a='pos'))
	# weight = scipy.linalg.solve(Sigma[:(N+2), :(N+2)], Sigma[N + 2][:(N+2)], assume_a='pos')[N + 1]

	# print(Sigma[N+1][N+2] - np.dot(Sigma[N+1][:(N+1)], scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+2][:(N+1)], assume_a='pos')))
	# print(scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+2][:(N+1)], assume_a='pos')[N] / precision)

	ans = Sigma[N+1][N+2] - np.dot(Sigma[N+1][:(N+1)], scipy.linalg.solve(Sigma[:(N+1), :(N+1)], Sigma[N+2][:(N+1)], assume_a='pos'))

	return ans


N = 2
D = 1

for i in range(1):

	X_initial = np.array([[2.0 - random.random()*4.0 for j in range(D)] for i in range(N)])
	err2 = np.zeros((N), dtype='double')
	for i in range(N):
		err2[i] = np.exp(3.0 - 6.0*random.random())
	X_target = np.array([0.0 for j in range(D)])

	X0 = np.array([2.0 - random.random()*4.0 for j in range(D)])
	X1 = np.array([2.0 - random.random()*4.0 for j in range(D)])
	precision = np.exp(3.0 - 6.0*random.random())

	res = np.zeros((2, 2), dtype='double')

	res[0][0] = calc(X_initial, err2, X0, precision, X0, X_target)
	res[0][1] = calc(X_initial, err2, X0, precision, X1, X_target)
	res[1][0] = calc(X_initial, err2, X1, precision, X0, X_target)
	res[1][1] = calc(X_initial, err2, X1, precision, X1, X_target)

	target0 = target(X_initial, err2, X0, precision, X_target)
	target1 = target(X_initial, err2, X1, precision, X_target)
	target2 = target(X_initial, err2, X_target, precision, X_target)

	print(res)
	# print(target(X_initial, err2, X0, precision, X_target))
	# print(target(X_initial, err2, X1, precision, X_target))

	if target1 > target0 and abs(res[0][0]) > abs(res[0][1]) and abs(res[1][1]) < abs(res[1][0]):
		print("CHECK 1 FAIL")
		print(res)
		print(target(X_initial, err2, X0, precision, X_target))
		print(target(X_initial, err2, X1, precision, X_target))
		break

	if target1 > target0 and target1 > target2 and abs(res[1][1]) < abs(res[1][0]):
		print("CHECK 2 FAIL")
		print(res)
		print(target0, target1, target2)
		break
