# This version is as follows:
# - It only supports Gaussian covariance function
# - The minimization procedure is a simple BFGS method.
# Input:
# - sequence l_0, l_1, ... , l_n
# - sequence r_k^2 or value r^2 or none
# - sequence p_k or value p
# - sequence w_k or none
# - matrices D and K of size s x s, where s = l_n + 1
# Key output:
# - new sequence l_0, l_1, ... , l_n
# - sequence r_k^2
# - sequence p_k
# - matrix X of size s x (n+1)
# Moreover, we implement an optimality tool that, for given index i, checks if there is room for improvement of evaluation.
# It extends 7th version.

# Sections:
# 0. Helpers
# 1. Main code
# 2. Optimality Checker
# 3. File reading and writing
# 4. Code execute

import os
import numpy as np
import scipy.linalg
from scipy.stats import ortho_group
from pathlib import Path

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.9f}\t'.format}, linewidth=200000, threshold=np.inf)

# Key variables
n = None
l = None
r2 = None
p = None
w = None
K = None
D = None

# Settings
MIN_CONTRIBUTION = 1.0/1e13
PRINTING = 0



# ===============================================================   HELPERS   ===============================================================

def softmax_reverse(u):
	global MIN_CONTRIBUTION
	sum = np.sum(u)

	v = []
	for x in u:
		if x > sum * MIN_CONTRIBUTION:
			v.append(np.log(x / sum))
		else:
			v.append(np.log(MIN_CONTRIBUTION) - 3.0)

	return np.array(v)


def collapse_onto_vector(X, precision):
	global n, l, r2, p

	v = np.array([], dtype='double')
	for k in range(n):
		for i in range(l[k], l[k+1]):
			v = np.concatenate((v, X[i][:(k+1)]))
	
	if r2 is None:
		v = np.concatenate((v, X[l[n]]))
	elif not hasattr(r2, "__len__"):
		v = np.concatenate((v, softmax_reverse(np.power(X[l[n]], 2))))

	if not hasattr(p, "__len__"):
		v = np.concatenate((v, softmax_reverse(precision)))
	else:
		for k in range(n):
			v = np.concatenate((v, softmax_reverse(precision[l[k]:l[k+1]])))

	return v


def softmax(v, scalar):
	global MIN_CONTRIBUTION

	u = []
	for y in v:
		if y < np.log(MIN_CONTRIBUTION):
			u.append(MIN_CONTRIBUTION)
		else:
			u.append(np.exp(y))

	u = np.array(u)
	u = u / np.sum(u) * scalar 
	return np.array(u)


def expand_vector(v):
	global n, l, r2, p

	X = np.zeros((l[n] + 1, n + 1), dtype='double')
	index = 0
	for k in range(n):
		for i in range(l[k], l[k+1]):
			X[i][:(k+1)] = v[index:(index+k+1)]
			index += k+1

	if r2 is None:
		X[l[n]] = v[index:(index+n+1)]
		index += n+1
	elif not hasattr(r2, "__len__"):
		X[l[n]] = np.sqrt(softmax(v[index:(index+n+1)], r2))
		index += n+1
	else:
		for k in range(n+1):
			X[l[n]][k] = np.sqrt(r2[k] - r2[k+1])

	precision = np.zeros((l[n]), dtype='double')
	if not hasattr(p, "__len__"):
		precision = softmax(v[index:index+l[n]], p)
		index += l[n]
	else:
		for k in range(n):
			precision[l[k]:l[k+1]] = softmax(v[index+l[k]:index+l[k+1]], p[k])
		index += l[n]

	return X, precision


# df/dv in terms of df/du, where u = softmax(v) * scalar
def softmax_gradient(u, u_gradient, scalar):
	sum = np.dot(u, u_gradient) / scalar
	return np.multiply(u, u_gradient) - sum * u


def gradient_to_v_gradient(X, precision, X_gradient, precision_gradient):
	global n, l, r2, p

	gradient = np.array([], dtype='double')
	for k in range(n):
		for i in range(l[k], l[k+1]):
			gradient = np.concatenate((gradient, X_gradient[i][:k+1]))

	if r2 is None:
		gradient = np.concatenate((gradient, X_gradient[l[n]]))
	elif not hasattr(r2, "__len__"):
		gradient = np.concatenate((gradient, softmax_gradient(np.power(X[l[n]], 2), np.divide(X_gradient[l[n]], 2.0 * X[l[n]]), r2)))

	if not hasattr(p, "__len__"):
		gradient = np.concatenate((gradient, softmax_gradient(precision, precision_gradient, p)))
	else:
		for k in range(n):
			gradient = np.concatenate((gradient, softmax_gradient(precision[l[k]:l[k+1]], precision_gradient[l[k]:l[k+1]], p[k])))

	return gradient


def minimize_wrapper(v):
	X, precision = expand_vector(v)
	ans, X_gradient, precision_gradient = calc(X, precision)
	return -ans, -gradient_to_v_gradient(X, precision, X_gradient, precision_gradient)



# ===============================================================   MAIN CODE   ===============================================================

# Calculates the result, and the derivatives with respect to all of X and precision.
def calc(X, precision):
	global PRINTING
	global n, l, r2, p, K, D

	# Calculate matrix K and ANS
	K_ = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	Sigma = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	for i0 in range(l[n] + 1):
		for i1 in range(l[n] + 1):
			Sigma[i0][i1] = np.exp(-(D[i0][i1] + np.dot(X[i0] - X[i1], X[i0] - X[i1])) / 2.0)
			K_[i0][i1] = Sigma[i0][i1] - np.exp(-(np.dot(X[i0], X[i0]) + np.dot(X[i1], X[i1])) / 2.0) * (np.exp(-D[i0][i1] / 2.0) - K[i0][i1])
	for i in range(l[n]):
		K_[i][i] += 1.0/precision[i]

	C = scipy.linalg.cholesky(K_, lower=True)


	# Calculate the answer and weights
	ans = 0.0
	weights = np.zeros((n+1), dtype='double')
	if hasattr(w, "__len__"):
		for k in range(n + 1):
			ans += w[k] * np.dot(C[l[n]][:l[k]], C[l[n]][:l[k]])
			weights[k] = w[k]
	else:
		for k in range(1, n+1):
			ans += X[l[n]][k] * np.sqrt(np.dot(C[l[n]][l[k-1]:l[k]], C[l[n]][l[k-1]:l[k]]))
			weights[k] = X[l[n]][k] / np.sqrt(np.dot(C[l[n]][l[k-1]:l[k]], C[l[n]][l[k-1]:l[k]]))
			if k != n:
				weights[k] -= X[l[n]][k+1] / np.sqrt(np.dot(C[l[n]][l[k]:l[k+1]], C[l[n]][l[k]:l[k+1]]))
		
	if PRINTING >= 1:
		print("ans:\t", ans) 
		print("weights:")
		print(weights)


	# Calculate conditioned matrices and other helpful matrices
	A_precision = np.zeros((n+1, l[n]), dtype='double')
	C_inverted = scipy.linalg.solve_triangular(C[:l[n], :l[n]], np.eye(l[n]), lower=True)
	P_inverted = np.diag(1.0 / precision[:l[n]])
	for k in range(n):
		A_precision[k+1] = A_precision[k] + C[l[n], l[k]:l[k+1]] @ C_inverted[l[k]:l[k+1], :l[n]]

	if PRINTING >= 3:
		print("A_precision:")
		print(A_precision)

	M = A_precision.T @ np.diag(weights) @ A_precision
	M_Sigma_X = np.multiply(M, Sigma[:l[n], :l[n]] + P_inverted) @ X[:l[n], :]
	u = np.multiply(weights @ A_precision, Sigma[l[n], :l[n]])


	# Calculate X_gradient and precision_gradient
	X_gradient = np.zeros((l[n]+1, n+1), dtype='double')
	for k in range(n):
		for i in range(l[n]):
			for j in range(k+1):
				X_gradient[i][j] = -M_Sigma_X[i][j] + u[i] * X[l[n]][j]

	if not hasattr(r2, "__len__"):
		sum = 0.0
		for k in range(1, n+1):
			sum += weights[k] * np.dot(C[l[n], :l[k]], C[l[n], :l[k]])
		for k in range(n+1):
			X_gradient[l[n]][k] = np.dot(u, X[:l[n], k]) - sum * X[l[n]][k]
			if not hasattr(w, "__len__") and k > 0:
				X_gradient[l[n]][k] += np.sqrt(np.dot(C[l[n], l[k-1]:l[k]], C[l[n], l[k-1]:l[k]]))


	precision_gradient = np.diag(M) / (2.0 * precision**2)
	if hasattr(w, "__len__"):
		X_gradient *= 2.0
		precision_gradient *= 2.0

	return ans, X_gradient, precision_gradient



# ==============================================================================================================================

n = 2
l = [0, 2, 4]
# w = None
w = [0.0, 0.1, 0.4]

precision0 = np.random.random((l[n]))
# p = np.sum(precision0)
p = np.zeros((n), dtype='double')
for k in range(n):
	p[k] = np.sum(precision0[l[k]:l[k+1]])



s = l[n]+1
d = 2
Z = np.random.normal(0.0, 1.0, (s, d))
D = np.zeros((s, s), dtype='double')
for i0 in range(s):
	for i1 in range(s):
		D[i0][i1] = np.dot(Z[i0] - Z[i1], Z[i0] - Z[i1])

K = np.diag(np.random.random((s)))
for i in range(s):
	K[i][i] = K[i][i] / (1.0 - K[i][i])
if s > 1:
	U = ortho_group.rvs(dim = s)
	K0 = U @ K @ U.T

X0 = np.zeros((l[n]+1, n+1), dtype='double')
for k in range(n):
	for i in range(l[k], l[k+1]):
		X0[i][:(k+1)] = np.random.normal(0.0, 1.0, (k+1))
X0[l[n]] = np.abs(np.random.normal(0.0, 1.0, (n+1)))
# r2 = np.dot(X0[l[n]], X0[l[n]])
for k in range(n+2):
	r2[k] = np.dot(X0[l[n]][k:], X0[l[n]][k:])



v0 = collapse_onto_vector(X0, precision0)

# (X1, precision1) = expand_vector(v0)
# print(precision0)
# print(precision1)
# print(X0)
# print(X1)

ans, v_gradient = minimize_wrapper(v0)
for i in range(len(v0)):
	v0[i] -= 0.001
	ANS0, _ = minimize_wrapper(v0)
	v0[i] += 0.002
	ANS1, _ = minimize_wrapper(v0)
	v0[i] -= 0.001

	# print(i, '\t', (ANS1 - ANS0)/0.002, '\t', v_gradient[i])
	print(i, '\t', ((ANS1 - ANS0)/0.002 - v_gradient[i]) / max(abs((ANS1 - ANS0)/0.002), abs(v_gradient[i])))



# ans, X_gradient, precision_gradient = calc(X0, precision)
# for k in range(n):
# 	for i in range(l[k], l[k+1]):
# 		for j in range(k+1):
# 			X0[i][j] -= 0.001
# 			ANS0, _, _ = calc(X0, precision)
# 			X0[i][j] += 0.002
# 			ANS1, _, _ = calc(X0, precision)
# 			X0[i][j] -= 0.001

# 			print(i, j, '\t', (ANS1 - ANS0)/0.002, '    \t', X_gradient[i][j])

# # PRINTING = 4
# for k in range(n+1):
# 	X0[l[n]][k] -= 0.001
# 	ANS0, _, _ = calc(X0, precision)
# 	X0[l[n]][k] += 0.002
# 	ANS1, _, _ = calc(X0, precision)
# 	X0[l[n]][k] -= 0.001

# 	print(l[n], k, '\t', (ANS1 - ANS0)/0.002, '    \t', X_gradient[l[n]][k])

# for i in range(l[n]):
# 	precision[i] -= 0.001
# 	ANS0, _, _ = calc(X0, precision)
# 	precision[i] += 0.002
# 	ANS1, _, _ = calc(X0, precision)
# 	precision[i] -= 0.001

# 	print(i, '\t', (ANS1 - ANS0)/0.002, '    \t', precision_gradient[i])
