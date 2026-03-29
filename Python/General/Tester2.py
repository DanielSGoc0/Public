# General Covariance
# Crude implementation, without derivatives
# This time we consider an added auxiliary (and fixed) distance matrix D.
import numpy as np
import scipy
from pathlib import Path
import scipy.linalg

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.9f}\t'.format}, linewidth=200000)

PRINTING = False
R2 = None
L = None
s = None
P = None
N = None
D = None
modifiable_S = None
modifiable = None

def get_v(X, S):
	global N, s, modifiable, modifiable_S

	v = []
	for i in range(N + 1):
		for j in range(s):
			if modifiable[i][j]:
				v.append(X[i][j])

	for i in range(len(modifiable_S)):
		if modifiable_S[i]:
			v.append(S[i])

	return v


def get_X(v, X0, S0):
	global N, s, modifiable, modifiable_S
	X = np.zeros((N + 1, s), dtype='double')

	k = 0
	for i in range(N + 1):
		for j in range(s):
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
	# return np.exp(-r2 / 2.0) + 1.0
	# return 4.0 * np.exp(-r2 * 4.0) + np.exp(-r2 / 4.0)
	return 16.0 * np.exp(-r2 * 8.0) + 1.0 * np.exp(-r2 / 16.0)


def calc(args):
	global PRINTING, R2, L, s, P, N, modifiable_S
	X = args[0]
	S = args[1]
	X = np.array(X)

	sigma2 = np.zeros((N + 1), dtype='double')

	for i in range(N):
		if not modifiable_S[i]:
			sigma2[i] = S[i]

	for k in range(s):
		W = 0.0
		for i in range(L[k], L[k + 1]):
			if modifiable_S[i]:
				W += np.exp(S[i])
		for i in range(L[k], L[k + 1]):
			if modifiable_S[i]:
				sigma2[i] = 1.0/(np.exp(S[i]) * P[k] / W) 



	X[N] *= np.sqrt(R2 / np.dot(X[N], X[N]))

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]) + D[i][j])
		Sigma[i][i] += sigma2[i]

	C = scipy.linalg.cholesky(Sigma, lower=True)

	ANS = 0.0
	for k in range(s):
		ANS += np.abs(X[N][k]) * np.sqrt(np.dot(C[N, L[k]:L[k + 1]], C[N, L[k]:L[k + 1]]))



	if PRINTING:
		print("ANS:", ANS)
		# h = np.sqrt(np.dot(X[N, L[0]:], X[N, L[0]:]) / (1.0 - np.dot(C[N, :L[0]], C[N, :L[0]])))
		# print("h:", h)

		f = np.zeros((s + 1), dtype='double')
		for k in range(s):
			f[k] = np.abs(X[N][k]) / np.sqrt(np.dot(C[N, L[k]:L[k + 1]], C[N, L[k]:L[k + 1]]))

		w = np.zeros((s), dtype='double')
		for k in range(s):
			w[k] = f[k] - f[k + 1]

		print("w:\t", w)
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

# We assume that L[0] = 0
def minimize():
	global PRINTING, R2, L, s, P, N, D, modifiable, modifiable_S
	PRINTING = False
	L = [0, 1, 2, 3]
	s = len(L) - 1
	N = L[-1]
	R2 = 1.0
	READ = True

	P = []
	for k in range(s):
		w = np.random.random()
		w = (w / (1.0 - w))**3
		P.append(w)



	INPUT = None
	X0 = np.zeros((N + 1, s), dtype='double')
	S0 = np.random.random((N))
	D = np.zeros((N + 1, N + 1), dtype='double')

	if READ:
		INPUT = read_from_file("in.txt")
		D = np.array(INPUT[0])
		X0 = np.array(INPUT[1])
		S0 = np.array(INPUT[2])
	else:
		for _ in range(1):
			x = np.random.normal(0.0, 1.0, (N + 1))
			for i in range(N + 1):
				for j in range(N + 1):
					D[i][j] += (x[i] - x[j])**2


	modifiable = [[False for j in range(s)] for i in range(N + 1)]

	for k in range(s):
		for i in range(L[k], L[k + 1]):
			for j in range(k):
				X0[i][j] = np.random.normal(0.0, 1.0, (1))
				modifiable[i][j] = True

	for j in range(s):
		X0[N][j] = np.random.normal(0.0, 1.0, (1))
		modifiable[N][j] = True


	X0[1][0] += 1.304687565
	X0[2][0] += 0.657444652
	X0[2][1] += 0.406976404


	modifiable_S = [False for i in range(N)]
	for i in range(0, N):
		if not READ:
			modifiable_S[i] = False
			w = np.random.random()
			w = (w / (1.0 - w))
			S0[i] = w



	fun = lambda v: -calc(get_X(v, X0, S0))

	RES = scipy.optimize.minimize(fun, get_v(X0, S0), options={"gtol": 0.0000000000000000000001})
	if PRINTING:
		print("=================================================================")
	# PRINTING = True
	ANS = fun(RES.x)

	X1 = get_X(RES.x, X0, S0)[0]
	X1[N] *= np.sqrt(R2 / np.dot(X1[N], X1[N]))
	# with open("out.txt", "a") as f:
	# 	f.write(str(np.array(D)))
	# 	f.write("\n")
	# 	f.write(str(X1))
	# 	f.write("\n")
	# 	f.write(str(get_X(RES.x, X0, S0)[1]))
	# 	f.write("\n\n")

	return (-ANS, X1)

BEST_X = None
BEST_ANS = 0.0

for k in range(100000):
	(ANS, X) = minimize()
	if ANS > BEST_ANS:
		BEST_ANS = ANS
		BEST_X = X

	if k % 100 == 0:
		print()
		print(k)
		print(BEST_ANS)
		print(BEST_X)
