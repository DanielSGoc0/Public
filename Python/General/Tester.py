# General Covariance
# Crude implementation, without derivatives
import numpy as np
import scipy
import random
from pathlib import Path

import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.9f}\t'.format}, linewidth=200000)

PRINTING = False
R2 = None
L = None
P = None
N = None
modifiable_S = None
modifiable = None

def get_v(X, S):
	global N, modifiable, modifiable_S

	v = []
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				v.append(X[i][j])

	for i in range(len(modifiable_S)):
		if modifiable_S[i]:
			v.append(S[i])

	return v


def get_X(v, X0, S0):
	global N, modifiable, modifiable_S
	X = np.zeros((N + 1, N + 1), dtype='double')

	k = 0
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				X[i][j] = v[k]
				k += 1
			else:
				X[i][j] = X0[i][j]
		
	S = np.zeros(N, dtype='double')
	for i in range(N):
		if modifiable_S[i]:
			S[i] = v[k]
			k += 1
		else:
			S[i] = S0[i]

	return (X, S)


# ===============================================================================

def phi(r2):
	# return np.exp(-r2 / 2.0)
	return np.exp(-r2 / 2.0) + 1.0
	# return 4.0 * np.exp(-r2 * 4.0) + np.exp(-r2 / 4.0)
	# return 16.0 * np.exp(-r2 * 8.0) + 1.0 * np.exp(-r2 / 16.0)

# def flipped(args):
# 	global PRINTING, R2, L, P, N, modifiable_S
# 	X = args[0]
# 	S = args[1]
# 	X = np.array(X)

# 	sigma2 = np.zeros((N + 1), dtype='double')

# 	for i in range(N):
# 		if not modifiable_S[i]:
# 			sigma2[i] = S[i]

# 	for k in range(len(L) - 1):
# 		W = 0.0
# 		for i in range(L[k], L[k + 1]):
# 			if modifiable_S[i]:
# 				W += np.exp(S[i])
# 		for i in range(L[k], L[k + 1]):
# 			if modifiable_S[i]:
# 				sigma2[i] = 1.0/(np.exp(S[i]) * P[k] / W) 


# 	# Now calculate the Sigma and C.
# 	Sigma = np.zeros((N + 1, N + 1), dtype='double')
# 	for i in range(N + 1):
# 		for j in range(N + 1):
# 			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))
# 		Sigma[i][i] += sigma2[i]

# 	C = scipy.linalg.cholesky(Sigma, lower=True)


# 	X2 = np.copy(X)
# 	for k in range(len(L) - 1):
# 	# for k in range(len(L) - 2, len(L) - 1):
# 		u = C[N, L[k]:L[k+1]] / np.sqrt(np.dot(C[N, L[k]:L[k+1]], C[N, L[k]:L[k+1]]))

# 		for i in range(L[k + 1], N + 1):
# 			# If inner product is negative, flip around the C vector.
# 			if np.dot(X2[i, L[k]:L[k+1]], u) < 0.0:
# 				s = 2.0 * np.dot(X2[i, L[k]:L[k+1]], u)
# 				X2[i, L[k]:L[k+1]] -= s * u

# 		# Now recalculate Sigma and C.
# 		for i in range(N + 1):
# 			for j in range(N + 1):
# 				Sigma[i][j] = phi(np.dot(X2[i] - X2[j], X2[i] - X2[j]))
# 			Sigma[i][i] += sigma2[i]

# 		C = scipy.linalg.cholesky(Sigma, lower=True)
# 		# print("next iter:")
# 		# print(C)
# 		# print()
				
# 	return X2



def calc(args):
	global PRINTING, R2, L, P, N, modifiable_S
	X = args[0]
	S = args[1]
	X = np.array(X)

	sigma2 = np.zeros((N + 1), dtype='double')

	for i in range(N):
		if not modifiable_S[i]:
			sigma2[i] = S[i]

	for k in range(len(L) - 1):
		W = 0.0
		for i in range(L[k], L[k + 1]):
			if modifiable_S[i]:
				W += np.exp(S[i])
		for i in range(L[k], L[k + 1]):
			if modifiable_S[i]:
				sigma2[i] = 1.0/(np.exp(S[i]) * P[k] / W) 



	X[N, L[0]:N] *= np.sqrt(R2 / np.dot(X[N, L[0]:N], X[N, L[0]:N]))

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))
		Sigma[i][i] += sigma2[i]

	C = scipy.linalg.cholesky(Sigma, lower=True)
	# print()
	# print(X)
	# print(C)
	
	ANS = 0.0
	for k in range(len(L) - 1):
		ANS += np.sqrt(np.dot(X[N, L[k]:L[k + 1]], X[N, L[k]:L[k + 1]]) * np.dot(C[N, L[k]:L[k + 1]], C[N, L[k]:L[k + 1]]))
	ANS = 0.0
	for k in range(len(L) - 1):
		if L[k + 1] - L[k] == 1:
			if X[N][L[k]] * C[N][L[k]] < 0.0:
				for i in range(N + 1):
					X[i][L[k]] *= -1
		else:
			u = X[N, L[k]:L[k + 1]] / np.sqrt(np.dot(X[N, L[k]:L[k + 1]], X[N, L[k]:L[k + 1]]))
			u -= C[N, L[k]:L[k + 1]] / np.sqrt(np.dot(C[N, L[k]:L[k + 1]], C[N, L[k]:L[k + 1]]))
			u /= np.sqrt(np.dot(u, u))
			
			for i in range(N + 1):
				X[N, L[k]:L[k + 1]] -= 2.0 * np.dot(X[N, L[k]:L[k + 1]], u) * u

		ANS += np.dot(X[N, L[k]:L[k + 1]], C[N, L[k]:L[k + 1]])
	# ANS = np.dot(C[N, L[0]:N], X[N, L[0]:N])


	if PRINTING:
		print("ANS:", ANS)
		# h = np.sqrt(np.dot(X[N, L[0]:], X[N, L[0]:]) / (1.0 - np.dot(C[N, :L[0]], C[N, :L[0]])))
		# print("h:", h)

		for k in range(N):
			Sigma[k][k] -= sigma2[k]
		# print(np.linalg.det(Sigma[(N-1):(N+1), (N-2):N]))


		f = np.zeros((len(L)), dtype='double')
		for k in range(len(L) - 1):
			f[k] = np.sqrt(np.dot(X[N, L[k]:L[k + 1]], X[N, L[k]:L[k + 1]]) / np.dot(C[N, L[k]:L[k + 1]], C[N, L[k]:L[k + 1]]))

		w = np.zeros((len(L) - 1), dtype='double')
		for k in range(len(L) - 1):
			w[k] = f[k] - f[k + 1]

		print("w:\t", np.array(w))
		# print("p\t:", P)
		print("e:\t", sigma2)
		print(X)

	return ANS


# =============================================================================================================================

def read_number(txt, cursor):
	x = 0
	negative = False
	after_decimal = -1
	while True:
		c = txt[cursor]
		if c == ' ' or c == '\t' or c == '\n' or c == ']':
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
		if c == ' ' or c == '\t' or c == '\n':
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


# =========================================================================================


# def test():
# 	global PRINTING, R2, L, P, N, modifiable, modifiable_S
# 	PRINTING = False
# 	L = [1, 2, 3]
# 	N = L[-1]
# 	R2 = 1.0

# 	P = []
# 	for k in range(len(L) - 1):
# 		s = np.random.random()
# 		s = (s / (1.0 - s))**3
# 		P.append(s)



# 	X0 = np.zeros((N + 1, N + 1), dtype='double')
# 	modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
	
# 	for k in range(len(L) - 1):
# 		for i in range(L[k], L[k + 1]):
# 			for j in range(L[k]):
# 				X0[i][j] = np.random.normal(0.0, 1.0, (1))
# 				if True:
# 				# if j >= L[0]:
# 					modifiable[i][j] = True

# 	for j in range(N):
# 		X0[N][j] = np.random.normal(0.0, 1.0, (1))
# 		# if True:
# 		if j >= L[0]:
# 			modifiable[N][j] = True
# 	X0[N, L[0]:N] *= np.sqrt(R2 / np.dot(X0[N, L[0]:N], X0[N, L[0]:N]))



# 	modifiable_S = [False for i in range(N)]

# 	S0 = np.random.random((N))
# 	for i in range(L[0]):
# 		S0[i] = 1000000000.0

# 	for i in range(L[0], N):
# 		modifiable_S[i] = True

# 	RES1 = calc((X0, S0))
	


# 	X1 = flipped((X0, S0))
# 	RES2 = calc((X1, S0))
# 	EPS = 0.00001
# 	return RES2 + EPS > RES1


def minimize():
	global PRINTING, R2, L, P, N, modifiable, modifiable_S
	PRINTING = True
	L = [2, 3, 5, 6]
	N = L[-1]
	R2 = 1.0
	READ = True

	P = []
	for k in range(len(L) - 1):
		s = np.random.random()
		s = (s / (1.0 - s))**3
		P.append(s)


	INPUT = None
	X0 = np.zeros((N + 1, N + 1), dtype='double')
	S0 = np.random.random((N))

	if READ:
		INPUT = read_from_file("in.txt")
		X0 = np.array(INPUT[0])
		S0 = np.array(INPUT[1])

	modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
	
	for i in range(L[0]):
		for j in range(L[0]):
			if not READ or j >= L[0]:
				X0[i][j] = np.random.normal(0.0, 1.0, (1))

	for k in range(len(L) - 1):
		for i in range(L[k], L[k + 1]):
			for j in range(L[k]):
				if not READ or j >= L[0]:
					X0[i][j] = np.random.normal(0.0, 1.0, (1))
				# if True:
				if j >= L[0]:
					modifiable[i][j] = True

	for j in range(N):
		if not READ or j >= L[0]:
			X0[N][j] = np.random.normal(0.0, 1.0, (1))
		# if True:
		if j >= L[0]:
			modifiable[N][j] = True



	modifiable_S = [False for i in range(N)]
	for i in range(L[0]):
		if not READ:
			s = np.random.random()
			s = (s / (1.0 - s))
			S0[i] = s
			# S0[i] = 1000000000.0



	# for i in range(L[0], N):
	# 	modifiable_S[i] = True
	for i in range(L[0], N):
		if not READ:
			modifiable_S[i] = False
			s = np.random.random()
			s = (s / (1.0 - s))
			S0[i] = s



	fun = lambda v: -calc(get_X(v, X0, S0))

	RES = scipy.optimize.minimize(fun, get_v(X0, S0), options={"gtol": 0.0000000000000000000001})
	print("=================================================================")
	# PRINTING = True
	fun(RES.x)

	with open("out.txt", "a") as f:
		f.write(str(get_X(RES.x, X0, S0)[0]))
		f.write("\n")
		f.write(str(get_X(RES.x, X0, S0)[1]))
		f.write("\n\n")

minimize()
# for k in range(1000):
# 	minimize()


# for k in range(10000):
# 	RES = test()
# 	if RES:
# 		print(k)
# 	else:
# 		print("damn")
# 		break
