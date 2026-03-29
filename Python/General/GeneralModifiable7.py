# This version is as follows:
# - It only supports fixed Gaussian covariance function
# - The problem is characterized by fixed initial data set (matrix K and D), fixed r^2 and fixed precisions per step. 
# - The minimization procedure is a simple BFGS method.
# It extends 6th version.

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
INF = 1e13

# ===============================================================================

def get_v(args):
	global n, l, modifiableX
	(X, R, P) = args

	v = []
	for i in range(l[n]):
		for j in range(n + 1):
			if modifiableX[i][j]:
				v.append(X[i][j])

	for r in R:
		v.append(r)
	for p in P:
		v.append(p)

	return v


def get_args(v):
	global n, l, modifiableX

	X = np.zeros((l[n] + 1, n + 1), dtype='double')
	index = 0
	for i in range(l[n]):
		for j in range(n + 1):
			if modifiableX[i][j]:
				X[i][j] = v[index]
				index += 1
			else:
				X[i][j] = 0.0

	R = np.zeros((n + 1), dtype='double')
	for j in range(n + 1):
		R[j] = v[index]
		index += 1

	P = np.zeros((l[n]), dtype='double')
	for i in range(l[n]):
		P[i] = v[index]
		index += 1

	return [X, R, P]


def reconstruct(args):
	global n, r2, l, p, modifiableX, K, D, INF
	X = args[0]
	R = args[1]
	P = args[2]

	# Rescale and obtain X_out and precisions p.
	s = 0.0
	max_R = np.max(R)

	for j in range(n + 1):
		if max_R - R[j] <= 30:
			s += np.exp(R[j] - max_R)
	s = np.sqrt(r2 / s)

	for j in range(n + 1):
		if max_R - R[j] <= 30:
			X[l[n]][j] = s * np.exp((R[j] - max_R) / 2.0)
		else:
			X[l[n]][j] = 0.0


	sigma2 = np.zeros((l[n]), dtype='double')
	if hasattr(p, "__len__"):
		for k in range(n):
			s = 0.0
			max_P = np.max(P[l[k]:l[k + 1]])

			for i in range(l[k], l[k + 1]):
				if max_P - P[i] <= 30:
					s += np.exp(P[i] - max_P)
			s = s / p[k]

			for i in range(l[k], l[k + 1]):
				if max_P - P[i] <= 30:
					sigma2[i] = s * np.exp(max_P - P[i])
				else:
					sigma2[i] = INF
	else:
		s = 0.0
		max_P = np.max(P)

		for i in range(l[n]):
			if max_P - P[i] <= 30:
				s += np.exp(P[i] - max_P)
		s = s / p

		for i in range(l[n]):
			if max_P - P[i] <= 30:
				sigma2[i] = s * np.exp(max_P - P[i])
			else:
				sigma2[i] = INF


	return (X, sigma2)


# ===============================================================================

def calc(X, sigma2, PRINTING):
	global n, r2, l, p, modifiableX, K, D

	# Calculate matrix K and ANS
	K_ = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	Sigma = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	for i0 in range(l[n] + 1):
		for i1 in range(l[n] + 1):
			Sigma[i0][i1] = np.exp(-(D[i0][i1] + np.dot(X[i0] - X[i1], X[i0] - X[i1])) / 2.0)
			K_[i0][i1] = Sigma[i0][i1] - np.exp(-(np.dot(X[i0], X[i0]) + np.dot(X[i1], X[i1])) / 2.0) * (np.exp(-D[i0][i1] / 2.0) - K[i0][i1])
	for i in range(l[n]):
		Sigma[i][i] += sigma2[i]
		K_[i][i] += sigma2[i]

	C = scipy.linalg.cholesky(K_, lower=True)

	ANS = 0.0
	for k in range(n):
		ANS += X[l[n]][k + 1] * np.sqrt(np.dot(C[l[n], l[k]:l[k + 1]], C[l[n], l[k]:l[k + 1]]))
	if PRINTING >= 1:
		print("ANS:")
		print(ANS)
 

	# ======   Additional variables   ======

	w = np.zeros((n + 1), dtype='double')
	for k in range(1, n + 1):
		w[k] = X[l[n]][k] / np.sqrt(np.dot(C[l[n], l[k - 1]:l[k]], C[l[n], l[k - 1]:l[k]]))
		if k < n:
			w[k] -= X[l[n]][k + 1] / np.sqrt(np.dot(C[l[n], l[k]:l[k + 1]], C[l[n], l[k]:l[k + 1]]))
	if PRINTING >= 2:
		print("w:")
		print(w)

	A = np.zeros((n + 1, l[n]), dtype='double') # A[k][i] non-zero if i < l[k], and has an additional factor of p[i]
	Cinv = scipy.linalg.solve_triangular(C[:l[n], :l[n]], np.eye(l[n]), lower=True)
	# for k in range(n + 1):
	# 	A[k, :l[k]] = C[l[n], :l[k]] @ Cinv[:l[k], :l[k]]
	for k in range(1, n + 1):
		A[k] = A[k - 1] + C[l[n], l[k-1]:l[k]] @ Cinv[l[k-1]:l[k], :l[n]]
	if PRINTING >= 3:
		print("A:")
		print(A)

	M = A.T @ np.diag(w) @ A
	M_Sigma_X = np.multiply(M, Sigma[:l[n], :l[n]]) @ X[:l[n], :]
	u = np.multiply(w @ A, Sigma[l[n], :l[n]])


	if PRINTING >= 4:
		print("r2:\t", r2)
		print("K:\t", K[l[n]][l[n]])
		print("X / XN:")
		I = []
		for i in range(l[n]):
			if sigma2[i] < INF/100.0:
				I.append(i)
		print(np.linalg.solve(np.multiply(M[np.ix_(I, I)], Sigma[np.ix_(I, I)]), u[I]))

	# ======   Calculating gradient   ======

	gradient = []
	for i in range(l[n]):
		for j in range(n + 1):
			if modifiableX[i][j]:
				gradient.append(-M_Sigma_X[i][j] + u[i] * X[l[n]][j])


	r2_gradient = np.zeros((n + 1), dtype='double')
	for j in range(n + 1):
		if j > 0:
			r2_gradient[j] = np.sqrt(np.dot(C[l[n], l[j - 1]:l[j]], C[l[n], l[j - 1]:l[j]]))
		r2_gradient[j] += np.dot(u, X[:l[n], j])

	for j0 in range(n + 1):
		partial = 1.0/2.0 * r2_gradient[j0] * X[l[n]][j0]
		for j1 in range(n + 1):
			partial -= 1.0/2.0 * r2_gradient[j1] * X[l[n]][j0]**2 * X[l[n]][j1] / r2
		gradient.append(partial)


	if hasattr(p, "__len__"):
		for k in range(n):
			for i0 in range(l[k], l[k + 1]):
				partial = M[i0][i0] * sigma2[i0] / 2.0
				for i1 in range(l[k], l[k + 1]):
					partial -= M[i1][i1] * sigma2[i1] / (2.0 * p[k] * sigma2[i0])
				gradient.append(partial)
	else:
		for i0 in range(l[n]):
			partial = M[i0][i0] * sigma2[i0] / 2.0
			for i1 in range(l[n]):
				partial -= M[i1][i1] * sigma2[i1] / (2.0 * p * sigma2[i0])
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
	np.random.seed(seed) # seed can be None

	n = 4
	# r2 = 0.4
	l = [0, 2, 4, 6, 8]
	# l = [0, 1, 2, 3]
	# p = [1.0]
	# n = 2
	# r2 = 1.0
	# l = [0, 1, 2]
	# p = [1.0, 1.0]
	# n = 3
	# r2 = 1.0
	# l = [0, 2, 4, 6]
	# p = [1.0, 1.0, 1.0]

	r2 = np.random.random()
	r2 = r2 / (1.0 - r2)
	w = np.random.random()
	w = w / (1.0 - w)
	# p = []
	# for k in range(n):
	# 	v = np.random.random()
	# 	v = v / (1.0 - v)
	# 	p.append(v * w)
	p = np.random.random()
	p = p / (1.0 - p)

	modifiableX = [[False for j in range(l[n] + 1)] for i in range(l[n] + 1)]
	for j in range(n + 1):
		for i in range(l[j], l[n]):
			modifiableX[i][j] = True

	D = np.zeros((l[n] + 1, l[n] + 1), dtype='double')
	K = np.zeros((l[n] + 1, l[n] + 1), dtype='double')

	d = 1
	X0 = np.random.normal(0.0, 0.3, (l[1], d))
	for i0 in range(l[1]):
		for i1 in range(l[1]):
			# D[i0][i1] = np.dot(X0[i0] - X0[i1], X0[i0] - X0[i1])
			D[i0][i1] = 0.0

	K0 = np.diag(np.random.random((l[1])))
	for i in range(l[1]):
		K0[i][i] = K0[i][i] / (1.0 - K0[i][i])
	if l[1] > 1:
		U = ortho_group.rvs(dim = l[1])
		K0 = U @ K0 @ U.T

	for i0 in range(l[1]):
		for i1 in range(l[1]):
			K0[i0][i1] = K0[0][0]

	# s = np.random.random()
	# for i in range(l[1]):
	# 	K0[i][0] *= s
	# 	K0[0][i] *= s

	for i0 in range(l[n] + 1):
		for i1 in range(l[n] + 1):
			D[i0][i1] = D[i0 % l[1]][i1 % l[1]]
			K[i0][i1] = K0[i0 % l[1]][i1 % l[1]]


def generate_from_file():
	global n, r2, l, p, modifiableX, K, D, INF
	info = read_from_file("in.txt")
	
	n = info[0]
	r2 = info[1]
	l = info[2]
	p = info[3]
	K = info[4]
	D = info[5]
	for k in range(n):
		sum = 0.0
		for i in range(l[k], l[k + 1]):
			sum += p[i]
		p[k] = sum

	modifiableX = [[False for j in range(l[n] + 1)] for i in range(l[n] + 1)]
	for j in range(n + 1):
		for i in range(l[j], l[n]):
			modifiableX[i][j] = True


def start_random(seed):
	global n, r2, l, p, modifiableX, K, D
	np.random.seed(seed) # seed can be None

	X_start = np.zeros((l[n] + 1, n + 1), dtype='double')
	for j in range(n + 1):
		for i in range(l[j], l[n]):
			X_start[i][j] = np.random.normal(0.0, 1.0)

	R_start = np.random.normal(0.0, 1.0, (n + 1))
	P_start = np.random.normal(0.0, 1.0, (l[n]))

	return (X_start, R_start, P_start)


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

