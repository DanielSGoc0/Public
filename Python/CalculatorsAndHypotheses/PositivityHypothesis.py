# Checking a hypothesis: is certain condition equivalent to being positive-definite?
# The answer is: YES
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import scipy.linalg

EPS = 0.00001

n = 4

def positive_definite(M):
	global EPS
	(evals, evecs) = np.linalg.eigh(M)
	POSITIVE = 1

	for v in evals:
		if v < -EPS:
			return -1
		elif v < EPS:
			POSITIVE = 0

	return POSITIVE

def hypothesis(M):
	global EPS, n
	e = np.array([1.0 - 2.0 * np.random.random() for j in range(n)])
	# e[0] = 1.0

	POSITIVE = 1

	w = 0.0
	for i in range(1, n + 1):
		v = e[:i] @ np.linalg.inv(M[:i, :i]) @ e[:i]
		# v = e[:i] @ M[:i, :i] @ e[:i]
		if v < w-EPS:
			return -1
		elif v < w+EPS:
			POSITIVE = 0
		w = v

	return POSITIVE



for i in range(1000000):
	A = np.array([[1.0 - 2.0 * np.random.random() for j in range(n)] for i in range(n)])
	A = A + A.T

	x = positive_definite(A)
	y = hypothesis(A)

	# if x == 1 and y == -1:
	if x * y == -1:
		print(A)
		print(x, y)
		break
