# This version is as follows:
# - It only supports fixed Gaussian covariance function
# - The problem is characterized by fixed initial data (matrix K and D), fixed delimiters l_k, fixed R^2 and fixed precisions. 
# - The minimization procedure is a simple BFGS method.

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
K = None
D = None


# ===============================================================================

def get_v(X):
	global n, l, modifiableX

	v = []
	for i in range(l[n] + 1):
		for j in range(n + 1):
			if modifiableX[i][j]:
				v.append(X[i][j])

	return v


def get_X(v):
	global n, l, modifiableX

	X = np.zeros((l[n] + 1, n + 1), dtype='double')
	index = 0
	for i in range(l[n] + 1):
		for j in range(n + 1):
			if modifiableX[i][j]:
				X[i][j] = v[index]
				index += 1
			else:
				X[i][j] = 0.0

	return X


# ===============================================================================

def calc(X, PRINTING):
	global n, r2, l, p, modifiableX, K, D
	X[l[n]][n] = 1.0
	scalar = np.sqrt(r2 / np.dot(X[l[n]], X[l[n]]))
	X[l[n]] *= scalar

	K_ = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	Sigma = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	for i0 in range(l[n] + 1):
		for i1 in range(l[n] + 1):
			Sigma[i0][i1] = np.exp(-(D[i0][i1] + np.dot(X[i0] - X[i1], X[i0] - X[i1])) / 2.0)
			K_[i0][i1] = Sigma[i0][i1] - np.exp(-(np.dot(X[i0], X[i0]) + np.dot(X[i1], X[i1])) / 2.0) * (np.exp(-D[i0][i1] / 2.0) - K[i0][i1])
	for i in range(l[n]):
		Sigma[i][i] += 1.0 / p[i]
		K_[i][i] += 1.0 / p[i]
	C = scipy.linalg.cholesky(K_, lower=True)

	ANS = 0.0
	for k in range(n):
		ANS += np.abs(X[l[n]][k + 1]) * np.sqrt(np.dot(C[l[n], l[k]:l[k + 1]], C[l[n], l[k]:l[k + 1]]))
	if PRINTING >= 1:
		print("ANS", ANS)
 

	# ======   Additional variables   ======

	w = np.zeros((n + 1), dtype='double')
	for k in range(1, n + 1):
		w[k] = np.abs(X[l[n]][k]) / np.sqrt(np.dot(C[l[n], l[k - 1]:l[k]], C[l[n], l[k - 1]:l[k]]))
		if k < n:
			w[k] -= np.abs(X[l[n]][k + 1]) / np.sqrt(np.dot(C[l[n], l[k]:l[k + 1]], C[l[n], l[k]:l[k + 1]]))
	if PRINTING >= 2:
		print(w)

	A = np.zeros((n + 1, l[n]), dtype='double')
	Cinv = scipy.linalg.solve_triangular(C[:l[n], :l[n]], np.eye(l[n]), lower=True)
	# for k in range(n + 1):
	# 	A[k, :l[k]] = C[l[n], :l[k]] @ Cinv[:l[k], :l[k]]
	for k in range(1, n + 1):
		A[k] = A[k - 1] + C[l[n], l[k-1]:l[k]] @ Cinv[l[k-1]:l[k], :l[n]]

	M = np.multiply(A.T @ np.diag(w) @ A, Sigma[:l[n], :l[n]])
	MX = M @ X[:l[n], :]
	u = np.multiply(w @ A, Sigma[l[n], :l[n]])


	# ======   Calculating gradient   ======

	last_gradient = np.zeros((n + 1), dtype='double')
	for j in range(n + 1):
		if j > 0:
			last_gradient[j] = np.sign(X[l[n]][j]) * np.sqrt(np.dot(C[l[n], l[j - 1]:l[j]], C[l[n], l[j - 1]:l[j]]))
		last_gradient[j] += np.dot(u, X[:l[n], j])

	gradient = []
	for i in range(l[n]):
		for j in range(n + 1):
			if modifiableX[i][j]:
				gradient.append(-MX[i][j] + u[i] * X[l[n]][j])

	for j in range(n + 1):
		if modifiableX[l[n]][j]:
			partial = scalar * last_gradient[j]
			for j1 in range(n + 1):
				partial -= scalar * X[l[n]][j] * X[l[n]][j1] / r2 * last_gradient[j1]
			gradient.append(partial)

	ANS = -ANS
	gradient = -np.array(gradient)
	return (ANS, np.array(gradient))



# ===============================================================================

def read_number(txt, cursor):
	x = 0
	negative = False
	after_decimal = -1
	while True:
		c = txt[cursor]
		if c == ' ' or c == '\t' or c == '\n' or c == ']' or c == ',':
			break
		elif c == '-':
			cursor += 1
			negative = True
		elif c == '.':
			cursor += 1
			after_decimal = 0
			continue
		else:
			cursor += 1
			x *= 10
			x += ord(c) - ord('0')
		if after_decimal >= 0:
			after_decimal += 1
	if after_decimal >= 0:
		x = float(x)
		x /= 10**after_decimal
	if negative:
		x = -x
	return (x, cursor)

def create_array(txt, cursor):
	RES = []
	while True:
		if cursor == len(txt):
			return (RES, cursor)
		c = txt[cursor]
		if c == ' ' or c == '\t' or c == '\n' or c == ',':
			cursor += 1
			continue
		elif c == '[':
			cursor += 1
			output = create_array(txt, cursor)
			RES.append(output[0])
			cursor = output[1]
		elif c == ']':
			cursor += 1
			return (RES, cursor)
		else:
			output = read_number(txt, cursor)
			RES.append(output[0])
			cursor = output[1]

def read_from_file(filename):
	txt = Path(filename).read_text()
	RES = create_array(txt, 0)[0]
	return RES


# ===============================================================================


def generate_random(seed):
	global n, r2, l, p, modifiableX, K, D

	n = 1
	r2 = 1.0
	l = [0, 1]
	p = [1.0]
	# n = 2
	# r2 = 1.0
	# l = [0, 1, 2]
	# p = [1.0, 1.0]
	# n = 3
	# r2 = 2.0
	# l = [0, 1, 2, 3]
	# p = [1.0, 1.0, 1.0]

	modifiableX = [[False for j in range(l[n] + 1)] for i in range(l[n] + 1)]
	for j in range(n + 1):
		for i in range(l[j], l[n] + 1):
			modifiableX[i][j] = True
	modifiableX[l[n]][n] = False


	np.random.seed(seed) # seed can be None

	D = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	K = np.zeros((l[n] + 1, l[n] + 1), dtype='double')

	d = 1
	X0 = np.random.normal(0.0, 1.0, (l[n] + 1, d))
	# X0 = np.random.uniform(0.0, 1.0, (l[n] + 1, d))
	for i0 in range(l[n] + 1):
		for i1 in range(l[n] + 1):
			D[i0][i1] = np.dot(X0[i0] - X0[i1], X0[i0] - X0[i1])
			# K[i0][i1] = np.exp(-D[i0][i1] / 2.0)
			# K[i0][i1] = 0.0

	K = np.diag(np.random.random((l[n] + 1)))
	for i in range(l[n] + 1):
		K[i][i] = K[i][i] / (1.0 - K[i][i])
	U = ortho_group.rvs(dim = l[n] + 1)
	K = U @ K @ U.T


def generate_from_file():
	global n, r2, l, p, modifiableX, K, D
	info = read_from_file("in.txt")
	
	n = info[0]
	r2 = info[1]
	l = info[2]
	p = info[3]
	K = info[4]
	D = info[5]

	modifiableX = [[False for j in range(l[n] + 1)] for i in range(l[n] + 1)]
	for j in range(n + 1):
		for i in range(l[j], l[n] + 1):
			modifiableX[i][j] = True
	modifiableX[l[n]][n] = False


def start_random(seed):
	global n, r2, l, p, modifiableX, K, D
	np.random.seed(seed) # seed can be None

	X_start = np.zeros((l[n] + 1, n + 1), dtype='double')
	for j in range(n + 1):
		for i in range(l[j], l[n] + 1):
			X_start[i][j] = np.random.normal(0.0, 1.0, (1))
	X_start[l[n]][n] = 1.0

	return X_start


def save_to_file(ans, X):
	global n, r2, l, p, modifiableX, K, D

	with open("out.txt", "a") as f:
		f.write(str(n))
		f.write("\n")
		f.write(str(r2))
		f.write("\n")
		f.write(str(l))
		f.write("\n")
		f.write(str(p))
		f.write("\n")
		f.write(str(K))
		f.write("\n")
		f.write(str(D))
		f.write("\n")
		f.write(str(ans))
		f.write("\n")
		f.write(str(X))
		f.write("\n\n\n")


def minimize(seed):
	global n, r2, l, p, modifiableX, K, D

	fun = lambda v: calc(get_X(v), 0)
	RES = scipy.optimize.minimize(fun, get_v(start_random(seed)), jac=True, options={"gtol": 0.0000000000000000000001}, method='BFGS')
	# RES = scipy.optimize.minimize(fun, get_v(start_random(seed)), jac=False, options={"gtol": 0.0000000000000000000001}, method='BFGS')

	X_end = get_X(RES.x)
	X_end[l[n]][n] = 1.0
	X_end[l[n]] *= np.sqrt(r2 / np.dot(X_end[l[n]], X_end[l[n]]))
	# print(X_end)
		
	return (-RES.fun, X_end)



if False:
	# generate_from_file()
	generate_random(0)
	(ans, X) = minimize(45)

	save_to_file(ans, X)

	print(X)
	X[l[n]] /= X[l[n]][n]
	calc(X, 2)

else:
	generate_random(3)

	max_ans = 0.0
	max_X = None
	max_seed = None

	for seed in range(100, 1000):
		(ans, X) = minimize(seed)
		print(seed, '\t', ans)
		
		if ans > max_ans:
			max_ans = ans
			max_X = X
			max_seed = seed

	print("\n")
	print(max_seed)
	print(max_ans)
	print(max_X)
	max_X[l[n]] /= max_X[l[n]][n]
	calc(max_X, 2)
	# print("\n")
	# for i in range(l[n] + 1):
	# 	for j in range(n + 1):
	# 		if modifiableX[i][j] and abs(max_X[l[n]][j]) > 0.0001:
	# 			max_X[i][j] /= max_X[l[n]][j]
	# print(max_X)
