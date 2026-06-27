# This is an experiment.
# The oracle always returns a vector in direction of target,
# but errored by amount roughly proportional to the distance.
import numpy as np
import scipy.linalg
from scipy.stats import ortho_group
from pathlib import Path

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=200000, threshold=np.inf)

n = None
r2 = None
l = None
p = None
modifiableX = None
INF = 1e13

# ===============================================================================

def get_v(args):
	global n, l, modifiableX
	(X) = args

	v = []
	for i in range(l[n] + 1):
		for j in range(l[n] + 1):
			if modifiableX[i][j]:
				v.append(X[i][j])

	return v


def get_args(v):
	global n, l, modifiableX

	X = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	index = 0
	for i in range(l[n] + 1):
		for j in range(l[n] + 1):
			if modifiableX[i][j]:
				X[i][j] = v[index]
				index += 1
			else:
				X[i][j] = 0.0

	return [X]


def reconstruct(args):
	global n, l, modifiableX
	X = args[0]

	return X


# ===============================================================================

# We assume that e2(r2) / r2 is increasing!!!
def e2(r2):
	return r2 * 10.0

def calc(X, PRINTING):
	global n, r2, l, modifiableX

	# Calculate sequence r2_[k] and e2_[k]
	r2_ = np.zeros((n + 1), dtype='double')
	e2_ = np.zeros((n + 1), dtype='double')
	r2_[0] = r2
	e2_[0] = e2(r2)

	for k in range(n + 1):
		X_sum = 0.0
		for j in range(k):
			r2_[k] += X[k][j]**2 * e2_[j]
		r2_[k] += X_sum**2 * r2

	if PRINTING >= 1:
		print("ANS:")
		print(r2_[n])
 
	return r2_[n]


# ===============================================================================


def generate_random(seed):
	global n, r2, l, modifiableX
	np.random.seed(seed) # seed can be None

	n = 4
	r2 = 0.4
	l = [0, 2, 4, 6, 8]

	modifiableX = [[False for j in range(l[n] + 1)] for i in range(l[n] + 1)]
	for j in range(n + 1):
		for i in range(l[j], l[n]):
			modifiableX[i][j] = True


def start_random(seed):
	global n, r2, l, modifiableX
	np.random.seed(seed) # seed can be None

	X_start = np.zeros((l[n] + 1, n + 1), dtype='double')
	for j in range(n + 1):
		for i in range(l[j], l[n]):
			X_start[i][j] = np.random.normal(0.0, 1.0)

	return (X_start)


def save_to_file(ans, X, sigma2):
	global n, r2, l, p, modifiableX, K, D

	with open("out.txt", "a") as f:
		f.write(str(n))
		f.write("\n")
		f.write(str(r2))
		f.write("\n")
		f.write(str(l))
		f.write("\n")
		f.write(str(np.reciprocal(sigma2)))
		f.write("\n")
		f.write(str(K))
		f.write("\n")
		f.write(str(D))
		f.write("\n")
		f.write(str(ans))
		f.write("\n")
		f.write(str(X))
		f.write("\n\n\n")


def minimize_func(v):
	(X, sigma2) = reconstruct(get_args(v))
	return calc(X, sigma2, 0)

def minimize(seed):
	global n, r2, l, p, modifiableX, K, D

	RES = scipy.optimize.minimize(minimize_func, get_v(start_random(seed)), jac=True, options={"gtol": 0.0000000000000000000001}, method='BFGS')
	# RES = scipy.optimize.minimize(fun, get_v(start_random(seed)), jac=False, options={"gtol": 0.0000000000000000000001}, method='BFGS')

	(X_end, sigma2_end) = reconstruct(get_args(RES.x))
	# print(X_end)
		
	return (-RES.fun, X_end, sigma2_end)



if False:
	generate_from_file()
	# generate_random(3)
	(ans, X, sigma2_end) = minimize(0)

	save_to_file(ans, X, sigma2_end)

	print("\n")
	print(ans)
	print(X)
	print(np.reciprocal(sigma2_end))
	print()
	calc(X, sigma2_end, 3)
	# print(X)
	# calc(X, sigma2_end, 2)

else:
	generate_random(646)
	# generate_from_file()

	max_ans = 0.0
	max_X = None
	max_sigma2 = None
	max_seed = None

	# epsilon = 0.00000000001
	# seed = 0
	# best_count = 0
	# while seed < 100 or best_count < 20: # at most 2^(-20) chance of failure, assuming uniform distribution
	# 	(ans, X, sigma2) = minimize(seed)
	# 	print(seed, '\t', ans)
		
	# 	if ans > max_ans:
	# 		if ans > max_ans + epsilon:
	# 			best_count = 1
	# 		max_ans = ans
	# 		max_X = X
	# 		max_sigma2 = sigma2
	# 		max_seed = seed
	# 	elif ans > max_ans - epsilon:
	# 		best_count += 1

	# 	seed += 1

	for seed in range(100):
		(ans, X, sigma2) = minimize(seed)
		print(seed, '\t', ans)
		
		if ans > max_ans:
			max_ans = ans
			max_X = X
			max_sigma2 = sigma2
			max_seed = seed

	print("\n")
	print(max_seed)
	print(max_ans)
	print(max_X)
	print(np.reciprocal(max_sigma2))
	print()
	calc(max_X, max_sigma2, 4)

	save_to_file(max_ans, max_X, max_sigma2)

