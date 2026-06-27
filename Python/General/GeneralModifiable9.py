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

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.13f}\t'.format}, linewidth=200000, threshold=np.inf)
np.seterr(all='raise')

# Key variables
n = None
l = None
r2 = None
p = None
w = None
K = None
D = None

# Settings
PRINTING = 0
MAX_OPTIMALITY_CHECKS=100

MIN_CONTRIBUTION = 1.0/1e13

MIN_PRECISION = 1.0/1e6
MAX_DISTANCE2 = 30.0
MIN_DISTANCE2 = 0.001

OPTIMALITY_TESTS_AMOUNT = 1000
OPTIMALITY_TESTS_VARIANCE = 1.0
OPTIMALITY_TESTS_FACTOR = 1.001

PRECISION_SCALING_FACTOR = 0.1


# ===========================================================================================================================================
# ===============================================================   HELPERS   ===============================================================
# ===========================================================================================================================================

def softmax_reverse(u):
	global MIN_CONTRIBUTION

	u = np.asarray(u, dtype="double")
	total = np.sum(u)

	v = np.full(
		u.shape,
		np.log(MIN_CONTRIBUTION) - 3.0,
		dtype="double",
	)

	active = u > total * MIN_CONTRIBUTION
	v[active] = np.log(u[active] / total)

	return v


def softmax(v, scalar):
	global MIN_CONTRIBUTION

	v = np.asarray(v, dtype="double")
	max_v = np.max(v)
	clipped_v = np.maximum(v - max_v, np.log(MIN_CONTRIBUTION))

	u = np.exp(clipped_v)
	u *= scalar / np.sum(u)

	return u


# df/dv in terms of df/du, where u = softmax(v) * scalar
def softmax_gradient(u, u_gradient, scalar):
	weighted_average = np.dot(u, u_gradient) / scalar
	return u * (u_gradient - weighted_average)


def collapse_onto_vector(X, precision, limiter):
	global n, r2, p

	out = limiter[n]
	parts = []

	# Pack the active lower-triangular parts of X.
	for k in range(n):
		parts.append(
			X[limiter[k]:limiter[k + 1], :k + 1].reshape(-1)
		)

	# Pack the final row of X.
	if r2 is None:
		parts.append(X[out])

	elif not hasattr(r2, "__len__"):
		parts.append(
			softmax_reverse(np.square(X[out]))
		)

	# Pack precision.
	if not hasattr(p, "__len__"):
		parts.append(
			softmax_reverse(precision)
		)

	else:
		for k in range(n):
			parts.append(
				softmax_reverse(
					precision[limiter[k]:limiter[k + 1]]
				)
			)

	return np.concatenate(parts)


def expand_vector(v, limiter):
	global n, r2, p

	out = limiter[n]
	X = np.zeros((out + 1, n + 1), dtype="double")

	index = 0

	# Unpack the active lower-triangular parts of X.
	for k in range(n):
		block_size = limiter[k + 1] - limiter[k]
		parameter_count = block_size * (k + 1)

		X[limiter[k]:limiter[k + 1], :k + 1] = (
			v[index:index + parameter_count]
			.reshape(block_size, k + 1)
		)

		index += parameter_count

	# Unpack the final row of X.
	if r2 is None:
		X[out] = v[index:index + n + 1]
		index += n + 1

	elif not hasattr(r2, "__len__"):
		X[out] = np.sqrt(
			softmax(
				v[index:index + n + 1],
				r2,
			)
		)

		index += n + 1

	else:
		X[out] = np.sqrt(
			r2[:-1] - r2[1:]
		)

	# Unpack precision.
	precision = np.zeros(out, dtype="double")

	if not hasattr(p, "__len__"):
		precision = softmax(
			v[index:index + out],
			p,
		)

		index += out

	else:
		for k in range(n):
			block_size = limiter[k + 1] - limiter[k]

			precision[limiter[k]:limiter[k + 1]] = softmax(
				v[index:index + block_size],
				p[k],
			)

			index += block_size

	return X, precision


def gradient_to_v_gradient(
	X,
	precision,
	X_gradient,
	precision_gradient,
	limiter
):
	global n, r2, p

	out = limiter[n]
	parts = []

	# Pack the active components of the X gradient.
	for k in range(n):
		parts.append(
			X_gradient[limiter[k]:limiter[k + 1], :k + 1].reshape(-1)
		)

	# Transform and pack the gradient of the final row of X.
	if r2 is None:
		parts.append(X_gradient[out])

	elif not hasattr(r2, "__len__"):
		parts.append(
			softmax_gradient(
				np.square(X[out]),
				X_gradient[out] / (2.0 * X[out]),
				r2,
			)
		)

	# Transform and pack the precision gradient.
	if not hasattr(p, "__len__"):
		parts.append(
			softmax_gradient(
				precision,
				precision_gradient,
				p,
			)
		)

	else:
		for k in range(n):
			parts.append(
				softmax_gradient(
					precision[limiter[k]:limiter[k + 1]],
					precision_gradient[limiter[k]:limiter[k + 1]],
					p[k],
				)
			)

	return np.concatenate(parts)


def minimize_wrapper(v, limiter, indices):
	X, precision = expand_vector(v, limiter)

	ans, _, _, _, X_gradient, precision_gradient = calc(X, precision, limiter, indices)

	return -ans, -gradient_to_v_gradient(X, precision, X_gradient, precision_gradient, limiter)


def test_derivative(v, limiter, indices):
	_, v_gradient = minimize_wrapper(v, limiter, indices)
	for i in range(len(v)):
		v[i] -= 0.001
		ans0, _ = minimize_wrapper(v, limiter, indices)
		v[i] += 0.002
		ans1, _ = minimize_wrapper(v, limiter, indices)
		v[i] -= 0.001

		print(i, '\t', (ans1 - ans0)/0.002, '\t', v_gradient[i])



# =============================================================================================================================================
# ===============================================================   MAIN CODE   ===============================================================
# =============================================================================================================================================

# Calculates the ans and derivatives with respect to X and precision.
# Additionally returns weights and Cholesky lower factor of covariance matrix.
def calc(X, precision, limiter, indices):
	global PRINTING
	global n, r2, p, K, D

	out = limiter[n]


	# -----------------------------------------------------------
	# Calculate Sigma and K_
	# -----------------------------------------------------------

	squared_norms = np.sum(np.square(X), axis=1)
	squared_distances = squared_norms[:, None] + squared_norms[None, :] - 2.0 * X @ X.T

	D_factor = np.exp(-0.5 * D[np.ix_(indices, indices)])
	Sigma = D_factor * np.exp(-0.5 * squared_distances)
	origin_factor = np.exp(-0.5 * (squared_norms[:, None] + squared_norms[None, :]))
	K_ = Sigma - origin_factor * (D_factor - K[np.ix_(indices, indices)])

	diagonal_indices = np.arange(out)
	K_[diagonal_indices, diagonal_indices] += 1.0 / precision

	C = scipy.linalg.cholesky(K_, lower=True, check_finite=False)


	# -----------------------------------------------------------
	# Calculate answer and weights
	# -----------------------------------------------------------

	gains = np.concatenate((np.zeros(1), np.cumsum(np.square(C[out, :out]))))[np.asarray(limiter)]

	weights = np.zeros(n + 1, dtype="double")
	block_norms = np.zeros(n + 1, dtype="double")

	if hasattr(w, "__len__"):
		weights[:] = w
		ans = np.dot(weights, gains)

	else:
		for k in range(1, n + 1):
			block_norms[k] = np.linalg.norm(C[out, limiter[k - 1]:limiter[k]])

		ans = np.dot(X[out, 1:], block_norms[1:])

		weights[1:] = X[out, 1:] / block_norms[1:]
		weights[1:n] -= X[out, 2:] / block_norms[2:]
		
	if PRINTING >= 5:
		print("ans:\t", ans)
		print("weights:")
		print(weights)


	# -----------------------------------------------------------
	# Calculate A_precision and other matrices
	# -----------------------------------------------------------
	
	right_hand_sides = np.zeros((n + 1, out), dtype="double")
	for k in range(1, n + 1):
		right_hand_sides[k, :limiter[k]] = C[out, :limiter[k]]

	A_precision = scipy.linalg.solve_triangular(C[:out, :out].T, right_hand_sides.T, lower=False, check_finite=False).T


	print()
	print(np.log(MIN_CONTRIBUTION))
	print(precision)
	print(np.log(precision))
	print(A_precision)
	print(weights)
	M = A_precision.T @ (weights[:, None] * A_precision)
	print(M)
	print(np.log(Sigma))
	M_Sigma_X = np.multiply(M, Sigma[:out, :out] + np.diag(1.0 / precision)) @ X[:out, :]
	
	u = (weights @ A_precision) * Sigma[out, :out]

	if PRINTING >= 7:
		print("A_precision:")
		print(A_precision)


	# -----------------------------------------------------------
	# Calculate gradient
	# -----------------------------------------------------------

	X_gradient = np.zeros((out+1, n+1), dtype='double')
	X_gradient[:out, :] = -M_Sigma_X + u[:, None] * X[out][None, :]

	if not hasattr(r2, "__len__"):
		weighted_gain = np.dot(weights, gains)

		X_gradient[out] = u @ X[:out] - weighted_gain * X[out]

		if not hasattr(w, "__len__"):
			X_gradient[out, 1:] += block_norms[1:]


	precision_gradient = np.diag(M) / (2.0 * np.square(precision))

	if hasattr(w, "__len__"):
		X_gradient *= 2.0
		precision_gradient *= 2.0

	return ans, C, weights, A_precision, X_gradient, precision_gradient



# ======================================================================================================================================================
# ===============================================================   OPTIMALITY CHECKER   ===============================================================
# ======================================================================================================================================================

def precision_gradient_at_x(X, C, weights, x, k, index, limiter, indices):
	global PRINTING
	global n, K, D

	out = limiter[n]

	squared_distances = np.sum(np.square(X - x), axis=1)
	squared_norms_X = np.sum(np.square(X), axis=1)
	squared_norm_x = np.dot(x, x)

	K_x_ = (
		np.exp(-0.5 * squared_distances - 0.5 * D[index][indices])
		- np.exp(-0.5 * (squared_norm_x + squared_norms_X))
		* (np.exp(-0.5 * D[index][indices]) - K[index][indices])
	)

	C_x = scipy.linalg.solve_triangular(C, K_x_, lower=True, check_finite=False)

	inner_products = np.concatenate((np.zeros(1), np.cumsum(C[out, :out] * C_x[:out])))[np.asarray(limiter[(k+1):])]

	res = np.dot(weights[(k+1):], np.square(K_x_[out] - inner_products))

	if hasattr(w, "__len__"):
		return res
	else:
		return 0.5 * res
	

def optimality_checker(X, precision, limiter, indices):
	global PRINTING
	global OPTIMALITY_TESTS_AMOUNT
	global OPTIMALITY_TESTS_VARIANCE
	global OPTIMALITY_TESTS_FACTOR
	global n

	out = limiter[n]

	ans, C, weights, _, _, precision_gradient = calc(X, precision, limiter, indices)

	if PRINTING >= 0:
		print("ans:", ans)

	if PRINTING >= 1:
		print("weights:")
		print(weights)

		for k in range(n):
			for i in range(limiter[k], limiter[k+1]):
				print("test", k, indices[i], precision_gradient[i], precision_gradient_at_x(X, C, weights, X[i], k, indices[i], limiter, indices))

	m = np.max(precision_gradient)
	for k in range(n):
		if hasattr(p, "__len__"):
			m = np.max(precision_gradient[limiter[k]:limiter[k+1]])

		for index in range(l[k], l[k+1]):
			for _ in range(OPTIMALITY_TESTS_AMOUNT):
				x = np.concatenate((X[out, :(k+1)] + np.random.normal(0.0, 1.0/(k+1), (k+1)), np.zeros((n-k))))

				x_precision_gradient = precision_gradient_at_x(X, C, weights, x, k, index, limiter, indices)

				# if k == 1:
				# 	print(k, i, x_precision_gradient, m, x)

				if x_precision_gradient > m * OPTIMALITY_TESTS_FACTOR:
					if PRINTING >= 1:
						print("found", k, index, x_precision_gradient, m)
						print(x)
					return k, index, x

	return -1, -1, None


def clean(X, precision, limiter, indices):
	global MIN_PRECISION, MAX_DISTANCE2, MIN_DISTANCE2
	global n, D

	out = limiter[n]

	new_X = []
	new_precision = []
	new_limiter = [0]
	new_indices = []

	# Process each interval [limiter[k], limiter[k+1])
	for k in range(n):
		saved_i = []

		for i in range(limiter[k], limiter[k + 1]):
			if precision[i] < MIN_PRECISION:
				continue
			if np.dot(X[i] - X[out], X[i] - X[out]) + D[indices[i]][indices[out]] > MAX_DISTANCE2:
				continue
			saved_i.append(i)

		# Build clusters by single-linkage: rows connected by eps_dist edges.
		unused = set(saved_i)

		while unused:
			root = unused.pop()
			cluster = [root]
			stack = [root]

			while stack:
				i = stack.pop()
				close = []

				for j in unused:
					if indices[j] == indices[root] and np.dot(X[i] - X[j], X[i] - X[j]) <= MIN_DISTANCE2:
						close.append(j)

				for j in close:
					unused.remove(j)
					stack.append(j)
					cluster.append(j)

			avg_row = np.average(X[cluster], axis=0, weights=precision[cluster])

			new_X.append(avg_row)
			new_precision.append(sum(precision[cluster]))
			new_indices.append(indices[root])

		new_limiter.append(len(new_precision))

	# Append the final X[out]
	new_X.append(X[out])
	new_indices.append(indices[out])

	return np.array(new_X), np.array(new_precision), new_limiter, new_indices


def total_minimizer(X, precision, limiter, indices):
	global PRINTING
	global MAX_OPTIMALITY_CHECKS
	global OPTIMALITY_TESTS_AMOUNT
	global OPTIMALITY_TESTS_VARIANCE
	global PRECISION_RESCALING_FACTOR
	global n

	for attempt in range(MAX_OPTIMALITY_CHECKS):
		v0 = collapse_onto_vector(X, precision, limiter)
		res = scipy.optimize.minimize(minimize_wrapper, v0, args=(limiter, indices), jac=True, options={"gtol": 0.0000000000000000000001}, method='BFGS')
		X, precision = expand_vector(res.x, limiter)

		X, precision, limiter, indices = clean(X, precision, limiter, indices)

		if PRINTING >= 1:
			print()
			print("--------------------------------------------------------")
			print()
			print("X:")
			print(X)
			print("precision:")
			print(precision)
			print("limiter:")
			print(limiter)
			print("indices:")
			print(indices)

		k, index, x = optimality_checker(X, precision, limiter, indices)

		if k > -1:
			X = np.insert(X, limiter[k+1], x, axis=0)

			if hasattr(p, "__len__"):
				precision_sum = np.sum(precision[limiter[k]:limiter[k+1]])
				precision[limiter[k]:limiter[k+1]] *= (1.0 - PRECISION_SCALING_FACTOR / n)
				precision = np.insert(precision, limiter[k+1], precision_sum * PRECISION_SCALING_FACTOR / n)
			else:
				precision_sum = np.sum(precision)
				precision *= (1.0 - PRECISION_SCALING_FACTOR / n)
				precision = np.insert(precision, limiter[k+1], precision_sum * PRECISION_SCALING_FACTOR / n)

			for k1 in range(k+1, n+1):
				limiter[k1] += 1

			indices.insert(limiter[k+1]-1, index)

		else:
			break

	if PRINTING >= 1:
		print()
		print("--------------------------------------------------------")
		print()

	return X, precision, limiter, indices



# ================================================================================================================================================
# ===============================================================   CODE EXECUTE   ===============================================================
# ================================================================================================================================================

def next_seed(increase):
	SEED_FILE = "seed.txt"

	seed = 838
	if os.path.exists(SEED_FILE):
		with open(SEED_FILE, "r") as f:
			seed = int(f.read().strip())

	# Increment seed
	if increase:
		seed = seed + 1

	# Save updated seed
	with open(SEED_FILE, "w") as f:
		f.write(str(seed))

	return seed


def print_input():
	global n, r2, p, w, l, K, D

	print("--------------------------  input  --------------------------")
	print("n:", n)
	print("r2:", r2)
	print("p:", p)
	print("w:", w)
	print("l:", l)
	print("K:")
	print(K)
	print("D:")
	print(D)


def print_output(X, precision, limiter, indices):
	global n, r2, p, l, K, D

	print("--------------------------  output  --------------------------")
	print("limiter:")
	print(limiter)
	print("precision:")
	print(precision)
	print("X:")
	print(X)
	print("D:")
	print(D[np.ix_(indices, indices)])
	print("K:")
	print(K[np.ix_(indices, indices)])


def generate_random_start():
	global n, r2, p, l, K, D

	out = l[n]

	X = np.zeros((out + 1, n + 1), dtype='double')

	if r2 is None:
		X[out] = np.random.random((n + 1)) / (n + 1)
		
	elif not hasattr(r2, "__len__"):
		X[out] = np.random.random((n + 1))
		X[out] *= np.sqrt(r2) / np.linalg.norm(X[out])

	else:
		X[out] = X[out] = np.sqrt(r2[:-1] - r2[1:])

	for k in range(n):
		for i in range(l[k], l[k+1]):
			X[i, :(k + 1)] = X[out, :(k + 1)]

	precision = np.random.random((out))
	if not hasattr(p, "__len__"):
		precision *= p / np.sum(precision)

	else:
		for k in range(n):
			precision[l[k]:l[k + 1]] *= p[k] / np.sum(precision[l[k]:l[k+1]])

	limiter = [v for v in l]
	indices = [v for v in range(out + 1)]

	return X, precision, limiter, indices



# ================================================================================================================================================
# ===============================================================   CODE EXECUTE   ===============================================================
# ================================================================================================================================================

def generate_next():
	global PRINTING
	global n, r2, p, l, w, K, D

	seed = next_seed(False)
	if PRINTING >= 0:
		print("seed:", seed)
	np.random.seed(seed) # seed can be None

	n = 3
	l = [0, 2, 4, 6]
	p = np.random.random(n)
	# p = np.random.random()
	p = p / (1.0 - p)
	w = np.random.random(n + 1)
	w[0] = 0.0
	r2 = np.random.random(n + 2)
	r2 = np.sum(r2) - np.cumsum(r2)

	# s = l[n] + 1
	s = l[n] + 1
	D = np.zeros((s, s), dtype='double')
	Z = np.random.normal(0.0, 1.0, (s, 3))
	for k in range(n):
		Z[2*k] = Z[s - 1]
	for i0 in range(s):
		for i1 in range(s):
			D[i0][i1] = np.dot(Z[i0] - Z[i1], Z[i0] - Z[i1])


	# K = np.diag(np.random.random((out + 1)))
	# for i in range(out + 1):
	# 	K[i][i] = K[i][i] / (1.0 - K[i][i])
	# if out + 1 > 1:
	# 	U = ortho_group.rvs(dim = out + 1)
	# 	K = U @ K @ U.T
	# K = np.ones((out + 1, out + 1), dtype='double')
	K = np.exp(-0.5 * D)

	# factor = np.random.random()
	# factor = factor / (1.0 - factor)
	# K *= factor





PRINTING = 0
MAX_DISTANCE2 = 50.0
OPTIMALITY_TESTS_AMOUNT = 100
MAX_OPTIMALITY_CHECKS = 100

MIN_CONTRIBUTION = 1.0/1e13

while True:
	generate_next()

	X0, precision0, limiter0, indices0 = generate_random_start()

	X1, precision1, limiter1, indices1 = total_minimizer(X0, precision0, limiter0, indices0)

	_, _, _, A_precision, _, _ = calc(X1, precision1, limiter1, indices1)

	works = True
	# for k in range(n):
	# 	if limiter1[k + 1] - limiter1[k] > 1:
	# 		works = False
	for k in range(n):
		works_single = False
		max_precision = 0.0
		for i in range(limiter1[k], limiter1[k + 1]):
			all_positive = True
			for j in range(k + 1, n + 1):
				if A_precision[j][i] < 0:
					all_positive = False
					break
			if max_precision < precision1[i]:
				max_precision = precision1[i]
				# works_single = all_positive and (2*precision1[i] > np.sum(precision1[limiter1[k]:limiter1[k + 1]]))
				works_single = all_positive
			# if A_precision[k + 1][i] > 0:
			# 	works_single = True
			# 	break
			# if all_positive:
			# 	works_single = True
			# 	break
		if not works_single:
			works = False
		# if not works_single and limiter1[k] < limiter1[k + 1]:
		# 	works = False

	if not works:
		print_input()
		print_output(X1, precision1, limiter1, indices1)
		print("A_precision:")
		print(A_precision)
		break
