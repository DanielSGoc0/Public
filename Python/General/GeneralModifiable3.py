# General covariance
# The third and hopefully last version!

import numpy as np
import scipy
import random
import scipy.linalg
from pathlib import Path
import itertools
import time

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f} \t'.format}, linewidth=200000, threshold=np.inf)

modifiableX = None
modifiableS = None
Z = None
P = 10.0
H = 0.5
R2 = None

FIXED_COORDS = False
FIXED_POINT = False



# ===================================================================

def get_v(N, X, S):
	global FIXED_COORDS, modifiable, modifiableS

	v = []
	if FIXED_COORDS:
		v = X[N][:N]
	else:
		for i in range(N + 1):
			for j in range(N + 1):
				if modifiable[i][j]:
					v.append(X[i][j])

	for i in range(len(modifiableS)):
		if modifiableS[i]:
			v.append(S[i])

	return v


def get_args(v, N, M, X0, S0):
	global FIXED_COORDS, modifiable, modifiableS

	X = np.zeros((N + 1, N + 1), dtype='double')
	k = 0
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				if FIXED_COORDS:
					X[i][j] = v[j]
				else:
					X[i][j] = v[k]
					k += 1
			else:
				X[i][j] = X0[i][j]
	if FIXED_COORDS:
		k = N

	S = np.zeros(len(modifiableS), dtype='double')
	for i in range(len(modifiableS)):
		if modifiableS[i]:
			S[i] = v[k]
			k += 1
		else:
			S[i] = S0[i]

	return [N, M, X, S]


def get_err2(N, S):
	global P, modifiableS

	err2 = np.zeros((N + 1), dtype='double')
	SUM = 0.0
	for i in range(N):
		if modifiableS[i]:
			SUM += np.exp(S[i])
	for i in range(N):
		if modifiableS[i]:
			err2[i] = 1.0/(np.exp(S[i]) * P / SUM)
		else:
			err2[i] = S[i]

	return err2


def get_S(N, err2):
	global P, modifiableS

	S = np.zeros((N), dtype='double')
	for i in range(N):
		if modifiableS[i]:
			S[i] = -np.log(P * err2[i])
		else:
			S[i] = err2[i]

	return S



# ===============================================================================

# calculates matrix inverse of a_{max(i, j)} * b_{min(i, j)}
def matrix_inverse(a, b, N):
	M = np.zeros((N, N), dtype='double')
	M[0][0] = b[1] / b[0] / (a[0] * b[1] - a[1] * b[0])
	M[N - 1][N - 1] = a[N - 2] / a[N - 1] / (a[N - 2] * b[N - 1] - a[N - 1] * b[N - 2])

	for k in range(0, N - 1):
		M[k][k + 1] = 1.0 / (a[k + 1] * b[k] - a[k] * b[k + 1])
		M[k + 1][k] = M[k][k + 1]
	
	for k in range(1, N - 1):
		M[k][k] = (a[k - 1] * b[k + 1] - a[k + 1] * b[k - 1]) / ((a[k] * b[k + 1] - a[k + 1] * b[k]) * (a[k - 1] * b[k] - a[k] * b[k - 1]))

	return M

# Checks if matrix M is totally positive.
def check_totally_positive(M, N):
	ERR = 0.000000001
	# loop over all submatrix sizes
	for k in range(1, N + 1):
		# choose row indices
		for rows in itertools.combinations(range(N), k):
			# choose column indices
			for cols in itertools.combinations(range(N), k):
				submatrix = M[np.ix_(rows, cols)]
				det = np.linalg.det(submatrix)

				if det <= -ERR:
					print(submatrix)
					return False
				elif det < 0.0:
					print(k, np.log(-det))
	return True

def count_sign_changes(M, N):
	RES = []
	for j in range(N):
		amt = 0
		for i in range(N - 1):
			if M[i][j] * M[i + 1][j] < 0.0:
				amt += 1
		RES.append(amt)
	return RES

def dodgson(A, N):
	EPS = 0.00000001
	if N == 1:
		return A[0][0] > -EPS
	
	B = np.zeros((N - 1, N - 1), dtype='double')
	for i in range(N - 1):
		for j in range(N - 1):
			B[i][j] = A[i][j] * A[i + 1][j + 1] - A[i][j + 1] * A[i + 1][j]
			if B[i][j] < -EPS:
				return False
	print(B * 1000000000)
	return dodgson(B, N - 1)


# ===============================================================================

# The original covariance function phi(r^2 / 2).
# By Bernstein's Theorem, it must be an integral of exponentials.
def phi(x):
	return np.exp(-x)

# The derivative of phi(r^2 / 2) = phi(x) with respect to x.
def dphi(x):
	return -np.exp(-x)

# Initial data is controlled via matrix Z and array modifiable.
def calc(args):
	global P, Z, H, modifiable, modifiableS, FIXED_POINT
	N = args[0]
	OFFSET = args[1]
	X = np.array(args[2])
	S = np.array(args[3])
	err2 = get_err2(N, S)

	D = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i + 1):
			D[i][j] = np.dot(X[i] - X[j], X[i] - X[j])
			D[j][i] = D[i][j]

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(D[i][j] / 2.0)
		Sigma[i][i] += err2[i]

	C = scipy.linalg.cholesky(Sigma, lower=True)
	EVAL = np.dot(np.transpose(Z), X[N])

	if FIXED_POINT:
		ANS = np.dot(EVAL, C[N]) - 1.0 / (2.0 * H) * np.dot(X[N], X[N]) - H / 2.0
	else:
		ANS = np.dot(EVAL, C[N])
	# ANS += 2.0 * np.exp(-np.dot(X[N], X[N]) / 2.0)
	# print(X @ np.diag(1.0 / C[N]))


	# ======   Calculating gradient   ======

	# EVALC = np.dot(Z, C[N])

	K = np.zeros((N + 1, N + 1), dtype='double')
	# Kprime = np.zeros((N + 1, N + 1), dtype='double')
	# Kbis = np.zeros((N + 1, N + 1), dtype='double')
	# Ktris = np.zeros((N + 1, N + 1), dtype='double')
	# L = np.zeros((N + 1, N + 1), dtype='double')
	# Lprime = np.zeros((N + 1, N + 1), dtype='double')
	# Lbis = np.zeros((N + 1, N + 1), dtype='double')
	# Ltris = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			K[i][j] = C[N][max(i, j)] * EVAL[min(i, j)]
			# Kprime[i][j] = C[N][min(i, j)] * EVAL[max(i, j)]
			# Kbis[i][j] = (ANS * EVAL[max(i, j)] - C[N][max(i, j)]) * X[N][min(i, j)]
			# Ktris[i][j] = (ANS * EVAL[min(i, j)] - C[N][min(i, j)]) * X[N][max(i, j)]
			# L[i][j] = EVALC[max(i, j)] * X[N][min(i, j)]
			# Lprime[i][j] = EVALC[min(i, j)] * X[N][max(i, j)]
			# Lbis[i][j] = (ANS * X[N][max(i, j)] - EVALC[max(i, j)]) * X[N][min(i, j)]
			# Ltris[i][j] = (ANS * X[N][min(i, j)] - EVALC[min(i, j)]) * X[N][max(i, j)]

	M = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(K), lower=False)
	M = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(M), lower=False)

	# Mprime = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(Kprime), lower=False)
	# Mprime = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(Mprime), lower=False)

	# Mbis = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(Kbis), lower=False)
	# Mbis = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(Mbis), lower=False)

	# Mtris = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(Ktris), lower=False)
	# Mtris = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(Mtris), lower=False)
	

	# ======   printing   ======

	# dSigma = np.zeros((N + 1, N + 1), dtype='double')
	# for i in range(N + 1):
	# 	for j in range(N + 1):
	# 		dSigma[i][j] = Sigma[i][j]
	# for i in range(OFFSET, N):
	# 	dSigma[i][i] = 1.0

	MSigma = np.zeros((N + 1, N + 1), dtype='double')
	# MdSigma = np.zeros((N + 1, N + 1), dtype='double')
	# MprimedSigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			MSigma[i][j] = Sigma[i][j] * M[i][j]
			# MdSigma[i][j] = dSigma[i][j] * M[i][j]
			# MprimedSigma[i][j] = dSigma[i][j] * Mprime[i][j]


	# Isigma2 = np.zeros((N + 1, N + 1), dtype='double')
	# Isigma = np.zeros((N + 1, N + 1), dtype='double')
	# for i in range(N + 1):
	# 	Isigma2[i][i] = err2[i]
	# 	Isigma[i][i] = np.sqrt(err2[i])

	# Xsigma = X @ Isigma

	# In = np.zeros((N + 1, N + 1), dtype='double')
	# Im = np.zeros((N + 1, N + 1), dtype='double')
	# for i in range(OFFSET, N):
	# 	In[i][i] = np.sqrt(-1.0 / (err2[i] * M[i][i]))
	# 	Im[i][i] = err2[i] * Mprime[i][i]

	# eNeN = np.zeros((N + 1, N + 1), dtype='double')
	# eNeN[N][N] = 1.0
	# XNXN = np.zeros((N + 1, N + 1), dtype='double')
	# CNXN = np.zeros((N + 1, N + 1), dtype='double')

	# for i in range(N + 1):
	# 	for j in range(N + 1):
	# 		CNXN[i][j] = C[N][i] * X[N][j]
	# 		XNXN[i][j] = X[N][i] * X[N][j]

	# print(MprimedSigma[OFFSET:N, OFFSET:N])
	# print(In[OFFSET:N, OFFSET:N] @ MprimedSigma[OFFSET:N, OFFSET:N] @ In[OFFSET:N, OFFSET:N])
	# W = In[OFFSET:N, OFFSET:N] @ MprimedSigma[OFFSET:N, OFFSET:N] @ In[OFFSET:N, OFFSET:N]
	# print(np.linalg.inv(W))
	# (evals, evecs) = np.linalg.eigh(Mprime)
	# print(evals)
	# print(evecs)

	# print(check_totally_positive(Sigma[OFFSET:N, OFFSET:N], N - OFFSET))
	# print(check_totally_positive(dSigma[OFFSET:N, OFFSET:N], N - OFFSET))


	# m = M[N - 1][N - 1] * err2[N - 1]**2
	# Z = (Isigma2 @ MdSigma @ Isigma2 / m)[OFFSET:N, OFFSET:N]
	# a = np.diagonal(Z, offset=1)
	# print(np.diagonal(Z, offset=1))
	# Z = Sigma[OFFSET:N, OFFSET:N]
	# b = np.diagonal(Z, offset=1)
	# # Z = np.divide(np.ones((N + 1, N + 1), dtype='double'), np.multiply(Sigma, Sigma))[OFFSET:N, OFFSET:N]
	# # Z = np.divide(np.ones((N + 1, N + 1), dtype='double'), Sigma)[OFFSET:N, OFFSET:N]
	# print(np.diagonal(Z, offset=1))
	# print("OK")
	# print(((1-a) - 2*(1-b))/np.multiply((1-b), (1-b)))

	# A = np.transpose(X) @ (-Im + MdSigma - ANS * eNeN) @ X + Lprime
	# A = (-Im + MdSigma - ANS * eNeN)
	# A = np.transpose(X[:, 0:N]) @ (-Im + MdSigma - ANS * eNeN)[:, :] @ X[:, 0:N] + Lprime[0:N, 0:N]
	# A = np.transpose(X[:, 0:N]) @ (-Im + MdSigma - ANS * eNeN)[:, :] @ X[:, 0:N] + Lprime[0:N, 0:N]
	# A = np.transpose(X[:, 0:N]) @ (-Im + eNv + np.transpose(eNv) - MprimedSigma - ANS * eNeN)[:, :] @ X[:, 0:N] + Lprime[0:N, 0:N]
	# A = np.transpose(X[:, OFFSET:N]) @ (-Im + eNv + np.transpose(eNv) - MprimedSigma - ANS * eNeN)[:, :] @ X[:, OFFSET:N] + Lprime[OFFSET:N, OFFSET:N]

	# A = ANS * XNXN - np.transpose(X) @ (MprimedSigma + Im) @ X - L
	# A = ANS * XNXN[OFFSET:N, OFFSET:N] - np.transpose(X[:, OFFSET:N]) @ (MprimedSigma + Im) @ X[:, OFFSET:N] - L[OFFSET:N, OFFSET:N]
	# A = ANS * XNXN[OFFSET:N, OFFSET:N] - L[OFFSET:N, OFFSET:N]
	# A = Proj @ (Lprime[OFFSET:N, OFFSET:N] - np.transpose(X[:, OFFSET:N]) @ Im @ X[:, OFFSET:N]) @ Proj
	# A = Lprime[OFFSET:N, OFFSET:N] - np.transpose(X[:, OFFSET:N]) @ Im @ X[:, OFFSET:N]
	# A = L[OFFSET:N, OFFSET:N] + np.transpose(X[:, OFFSET:N]) @ Im @ X[:, OFFSET:N]
	# A = Lprime[OFFSET:N, OFFSET:N] + L[OFFSET:N, OFFSET:N]

	# A = Kprime[OFFSET:, OFFSET:]
	# A = A[OFFSET:, OFFSET:]

	# A = -np.transpose(X) @ MSigma @ X + Lprime - ANS * XNXN
	# A = np.transpose(X) @ MSigma @ X - Ltris
	# print(A)

	# (evals, evecs) = np.linalg.eigh(A)
	# print(evals)

	# for i in range(OFFSET - 1, N):
	# 	print(np.dot(MSigma[(i + 1):N, N], np.linalg.solve(MSigma[(i + 1):N, (i + 1):N], MSigma[(i + 1):N, N])), ANS - C[N][i]/X[N][i])


	# iter = lambda S: X[OFFSET:, OFFSET:].T @ np.multiply(Sigma[OFFSET:, OFFSET:], scipy.linalg.solve_triangular(C[OFFSET:, OFFSET:].T, scipy.linalg.solve_triangular(C[OFFSET:, OFFSET:].T, S, lower=False).T, lower=False)) @ X[OFFSET:, OFFSET:] + CNXN[OFFSET:, OFFSET:] + CNXN[OFFSET:, OFFSET:].T - ANS * XNXN[OFFSET:, OFFSET:] 

	# S = np.eye(N + 1 - OFFSET, dtype='double')
	# for k in range(10):
	# 	print(S)
	# 	S = iter(S)

	# for j in range(OFFSET, N):
	# 	print()
	# 	print(np.linalg.solve((Im + MprimedSigma)[j:N, j:N], w[j:N]))
	# 	# if j > OFFSET:
	# 	# 	print(X[j:N, j - 1] / X[N][j - 1])


	# for J in range(OFFSET + 1, N):
	# 	print("J =", J)

	# 	A = MprimedSigma[(J + 1):N, (J + 1):N]

	# 	for i in range(J + 1, N):
	# 		A[i - (J + 1)][i - (J + 1)] -= M[i][i] * err2[i]
		


	# 	A = In[(J + 1):N, (J + 1):N] @ MprimedSigma[(J + 1):N, (J + 1):N] @ In[(J + 1):N, (J + 1):N]
	# 	# print(np.linalg.inv(A))
	# 	(evals, evecs) = np.linalg.eigh(A)
	# 	print(evals)
	# 	# print(evecs[-1])

	# Tq = np.zeros((N - OFFSET - 1, N - OFFSET - 1), dtype='double')
	# Tr = np.zeros((N - OFFSET - 1, N - OFFSET - 1), dtype='double')
	# Jq = np.zeros((N - OFFSET, N - OFFSET - 1), dtype='double')
	# Jr = np.zeros((N - OFFSET, N - OFFSET - 1), dtype='double')
	# for i in range(N - OFFSET - 1):
	# 	Tq[i][i] = 1.0/(C[N][i + OFFSET + 1]/X[N][i + OFFSET + 1] - C[N][i + OFFSET]/X[N][i + OFFSET])
	# 	Tr[i][i] = 1.0/(X[N][i + OFFSET]/C[N][i + OFFSET] - X[N][i + OFFSET + 1]/C[N][i + OFFSET + 1])
	# 	Jq[i][i] = -1.0/X[N][i + OFFSET]
	# 	Jq[i + 1][i] = 1.0/X[N][i + OFFSET + 1]
	# 	Jr[i][i] = -1.0/C[N][i + OFFSET]
	# 	Jr[i + 1][i] = 1.0/C[N][i + OFFSET + 1]


	# inv = matrix_inverse(C[N][OFFSET:], X[N][OFFSET:], N + 1 - OFFSET)
	# inv = matrix_inverse(C[N][OFFSET:N], X[N][OFFSET:N], N - OFFSET)
	# inv = matrix_inverse(X[N][OFFSET:N], (ANS * X[N][OFFSET:N] - C[N][OFFSET:N]), N - OFFSET)

	# A = C[OFFSET:, OFFSET:] @ inv @ C[OFFSET:, OFFSET:].T
	# print(A)

	# print(np.linalg.inv(M[(OFFSET + 1):, (OFFSET + 1):]))

	# A = -Jq @ Tq @ Jq.T
	# print(A)
	# print(Jr @ Tr @ Jr.T)
	# print(inv)
	# print(np.linalg.inv(K[(OFFSET + 1):, (OFFSET + 1):]))
	# print(np.linalg.inv(K[OFFSET:N, OFFSET:N]) + Jq @ Tq @ Jq.T)
	# print(1.0/(C[N][N] * X[N][N - 1]))
	# print(-C[N][N-1]/(C[N][N]**2 * X[N][N-1]))

	# W = np.zeros((N - OFFSET + 1, N - OFFSET + 1), dtype='double')
	# W[:(N - OFFSET), :(N - OFFSET)] = -Jq @ Tq @ Jq.T
	# W[N - OFFSET - 1][N - OFFSET] = 1.0/(C[N][N] * X[N][N - 1])
	# W[N - OFFSET][N - OFFSET - 1] = 1.0/(C[N][N] * X[N][N - 1])
	# W[N - OFFSET][N - OFFSET] = -C[N][N-1]/(C[N][N]**2 * X[N][N-1])

	# print(np.linalg.inv(W) - K[OFFSET:, OFFSET:])

	# print(np.linalg.inv(K[(OFFSET + 1):, (OFFSET + 1):])[:(N - OFFSET - 1), :(N - OFFSET - 1)] + (Jq @ Tq @ Jq.T)[1:, 1:])
	# print(1.0/X[N][OFFSET + 1]**2 * 1.0/(C[N][OFFSET + 2]/X[N][OFFSET + 2] - C[N][OFFSET + 1]/X[N][OFFSET + 1]))
	# print(1.0/X[N][OFFSET + 1]**2 * 1.0/(C[N][OFFSET + 1]/X[N][OFFSET + 1] - C[N][OFFSET]/X[N][OFFSET]))

	# print(np.linalg.inv(Ltris[OFFSET:N, OFFSET:N]))
	# print(1/X[N][OFFSET]**2 * 1.0/(ANS - C[N][OFFSET]/X[N][OFFSET]))

	# Lp = np.zeros((N - OFFSET - 1), dtype='double')
	# Lq = np.zeros((N - OFFSET - 2), dtype='double')
	# OPERATOR = np.zeros((N - OFFSET, N - OFFSET), dtype='double')
	# for k in range(OFFSET, N - 1):
	# 	Lp[k - OFFSET] = 1.0/(C[N][k + 1] * X[N][k] - C[N][k] * X[N][k + 1])
	# 	OPERATOR[k - OFFSET][k - OFFSET] -= Lp[k - OFFSET]
	# 	OPERATOR[k - OFFSET + 1][k - OFFSET + 1] -= Lp[k - OFFSET]
	# 	OPERATOR[k - OFFSET + 1][k - OFFSET] += Lp[k - OFFSET]
	# 	OPERATOR[k - OFFSET][k - OFFSET + 1] += Lp[k - OFFSET]
	# for k in range(OFFSET + 1, N - 1):
	# 	Lq[k - OFFSET - 1] = (C[N][k + 1] - C[N][k]) * (X[N][k] - X[N][k - 1]) - (C[N][k] - C[N][k - 1]) * (X[N][k + 1] - X[N][k])
	# 	Lq[k - OFFSET - 1] /= (C[N][k + 1] * X[N][k] - C[N][k] * X[N][k + 1])
	# 	Lq[k - OFFSET - 1] /= (C[N][k] * X[N][k - 1] - C[N][k - 1] * X[N][k])
	# 	OPERATOR[k - OFFSET][k - OFFSET] += Lq[k - OFFSET - 1]

	# print(OPERATOR)

	# print(Lp[1]*(1.0 - X[N][OFFSET + 2]/X[N][OFFSET + 1]) - Lp[1])

	# C[N][k + 1] / X[N][k + 1] - C[N][k] / X[N][k] = c / (X[N][k + 1] * X[N][k])
	# w = 0.0
	# for k in range(OFFSET, N - 1):
	# 	w += 1.0 / (X[N][k + 1] * X[N][k])
	# 	# print(1.0/(C[N][k + 1] * X[N][k] - C[N][k] * X[N][k + 1]))
	# 	# print((C[N][k + 1] / X[N][k + 1] - C[N][k - 1] / X[N][k - 1]) * X[N][k]**2)
	# 	v = C[N][k + 1] / X[N][k + 1] - C[N][OFFSET] / X[N][OFFSET]
	# 	v = w/v
	# 	# v = C[N][k + 1] * X[N][k] - C[N][k] * X[N][k + 1]
	# 	# v = (C[N][k + 1] - C[N][k]) * (X[N][k - 1] - X[N][k]) - (C[N][k - 1] - C[N][k]) * (X[N][k + 1] - X[N][k])
	# 	# v /= C[N][k + 1] * X[N][k] - C[N][k] * X[N][k + 1]
	# 	# v /= C[N][k] * X[N][k - 1] - C[N][k - 1] * X[N][k]
	# 	# v = (C[N][k + 1] * X[N][k] - C[N][k] * X[N][k + 1])/(C[N][k] * X[N][k - 1] - C[N][k - 1] * X[N][k])
	# 	print(v)
	# 	# print(X[N][k] / err2[k])
	# 	# v -= 1.0
	# 	# v *= err2[k]**2
	# 	# print(v)
	# 	# print(MSigma[N][k + 1] * X[N][k + 1] * X[k + 1][k])
	# 	# print(MSigma[N, (k+2):N] @ (X[N][k + 1] * X[(k+2):N, k] - X[N][k] * X[(k+2):N, k+1]))
	# 	# print(X[k + 1][k] / np.sqrt(err2[k]))
	# 	# print(np.dot(X[N][k:], X[N][k:]) / np.log((X[N][k] * X[N][k - 1])/(C[N][k + 1] * C[N][k])))
	# print(C[N][OFFSET] / X[N][OFFSET], C[N][N - 1] / X[N][N - 1], w)

	

	# =========================================================

	X_sums = MSigma @ X

	gradient = []
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable[i][j]:
				partial = X_sums[i][j]

				if i == N:
					partial -= ANS * X[N][j]
					partial += np.dot(Z[j], C[N])
					if FIXED_POINT:
						partial -= X[N][j] / H
					# When we add summand c * e^(-||X_N||^2/2) or similar,
					# the contribution to gradient is accounted for in ANS.

				gradient.append(partial)

	
	for i in range(N):
		if modifiableS[i]:
			partial = -M[i][i] * err2[i] / 2.0

			for j in range(N):
				if modifiableS[j]:
					partial += M[j][j] * np.exp(S[i] - S[j]) / (2.0 * P)

			gradient.append(partial)

	print("ANS:", ANS)
	h = np.sqrt(np.dot(X[N, OFFSET:], X[N, OFFSET:]) / (1.0 - np.dot(C[N, :OFFSET], C[N, :OFFSET])))
	print("h:", h)
	w = [X[N][k] / C[N][k] - X[N][k + 1] / C[N][k + 1] for k in range(OFFSET, N)]
	print("w:", np.array(w))
	v = [X[N][k] * np.sqrt(err2[k]) - X[N][k + 1] * np.sqrt(err2[k + 1]) for k in range(OFFSET, N)]
	print("v:", np.array(v))
	print("X:", X[N, OFFSET:N])
	print("C:", C[N, OFFSET:N])
	print("p:", 1.0 / err2[OFFSET:N])
	m = np.zeros((N - OFFSET), dtype='double')
	for k in range(OFFSET, N):
		Sigma_k = np.copy(Sigma[k, :])
		Sigma_k[k] = 1.0
		for i in range(k, N):
			m[k - OFFSET] += w[i - OFFSET] * (Sigma[k][N] -  Sigma_k[:(i + 1)] @ np.linalg.inv(Sigma[:(i + 1), :(i + 1)]) @ Sigma[:(i + 1), N])**2
	print("m:", m)


	ANS = -ANS
	gradient = -np.array(gradient)
	return (ANS, np.array(gradient))



# =========================================================================================

# Returns True if final RES is non-negative.
def extend_arrays(N, n, s, X1, Sigma1, C1, err1, A, quotients, FULL):
	global R2
	RES = R2

	# Now the main loop
	for k in range(n, N):

		# First update Sigma1[k][:k] and C1[k][:(k + 1)]
		for j in range(k):
			Sigma1[k][j] = np.exp(-np.dot(X1[k] - X1[j], X1[k] - X1[j]) / 2.0)
			Sigma1[j][k] = Sigma1[k][j]
			C1[k][j] = (Sigma1[k][j] - np.dot(C1[k, :j], C1[j, :j])) / C1[j][j]

		# Next update Sigma1[k][k] and C1[k][k] and Sigma1[N][k] and C1[N][k]
		Sigma1[k][k] = 1.0 + err1[k]
		C1[k][k] = np.sqrt(max(Sigma1[k][k] - np.dot(C1[k, :k], C1[k, :k]), err1[k]))
		Sigma1[N][k] = np.exp(-(np.dot(X1[N, :k] - X1[k, :k], X1[N, :k] - X1[k, :k]) + RES) / 2.0)
		C1[N][k] = (Sigma1[N][k] - np.dot(C1[N, :k], C1[k, :k])) / C1[k][k]

		# Now that we have C1[N][k], we can update X1[N][k] and R2_res
		X1[N][k] = C1[N][k] / (quotients[k] + s)
		RES -= X1[N][k]**2

		# If X[N, n:N] already has modulus greater than R, then we stop (unless FULL is true)
		if (not FULL) and RES < 0.0:
			return (False, 0.0)

		# Finally, we update X1[(k + 1):N, k].
		X1[(k + 1):N, k] = -np.linalg.solve(A[(k + 1):N, (k + 1):N], A[(k + 1):N, N] * X1[N][k])

	return (True, RES)


# We assume X is fixed along X[0:(N + 1), 0:n].
# Furthermore P and R2 are fixed and given.
def iterate(N, n, X0, err0):
	global P, Z, R2
	X1 = np.copy(X0)

	Sigma0 = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma0[i][j] = np.exp(-np.dot(X0[i] - X0[j], X0[i] - X0[j]) / 2.0)
		Sigma0[i][i] += err0[i]
	C0 = scipy.linalg.cholesky(Sigma0, lower=True)

	EVALX = np.dot(np.transpose(Z), X0[N])

	K = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			K[i][j] = C0[N][max(i, j)] * EVALX[min(i, j)]
	M = scipy.linalg.solve_triangular(np.transpose(C0), np.transpose(K), lower=False)
	M = scipy.linalg.solve_triangular(np.transpose(C0), np.transpose(M), lower=False)

	A = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			A[i][j] = M[i][j] * Sigma0[i][j]

	# ==========================================================
	EPS = 10.0**(-10)

	# With A at hand, we can go on to derive the equations for X1.
	# First we extract new err. By assumption M is negative definite, so Mii <= 0
	err1 = np.copy(err0)

	# ====================   Updates of sigma2 are disabled!!!!   ====================
	# sum = 0.0
	# for k in range(n, N):
	# 	sum += np.sqrt(max(EPS, -M[k][k]))
	# sum /= P
	
	# for k in range(n, N):
	# 	err1[k] = sum / np.sqrt(max(EPS, -M[k][k]))


	# quotients is an increasing non-negative sequence, such that quotients[N-1] = F0.
	F0 = np.dot(C0[N], X0[N])
	quotients = np.zeros((N + 1), dtype='double')
	for k in range(n, N):
		quotients[k] = F0 + A[N, (k+1):N] @ np.linalg.inv(A[(k+1):N, (k+1):N]) @ A[(k+1):N, N]
	# print(quotients)

	# Now initialize Sigma1 and C1, on indices [0:n, 0:n] and [N, 0:n]
	Sigma1 = np.zeros((N + 1, N + 1), dtype='double')
	C1 = np.zeros((N + 1, N + 1), dtype='double')
	for j in range(n):
		for i in range(n):
			Sigma1[i][j] = Sigma0[i][j]
			C1[i][j] = C0[i][j]
		Sigma1[N][j] = Sigma0[N][j]
		C1[N][j] = C0[N][j]

	# Next update indices X1[n:N, 0:n]
	for j in range(n):
		X1[n:N, j] = -np.linalg.solve(A[n:N, n:N], A[n:N, :n] @ X1[:n, j] + A[n:N, N] * X1[N][j])
		# X1[n:N, j] = X0[n:N, j]


	# Time for the main binsearch
	min_s = -quotients[n] # At min_s, the modulus of X[N, n:N] is greater than R.
	max_s = 1.0 / np.sqrt(R2) # At max_s, the modulus of X[N, n:N] is smaller than R.

	# for s in np.linspace(min_s, max_s, 10):
	# 	(success, val) = extend_arrays(N, n, s, X1, Sigma1, C1, err1, A, quotients, False)
	# 	print("s =", s, "   success =", success)

	right_sequence = []

	(success, val) = extend_arrays(N, n, max_s, X1, Sigma1, C1, err1, A, quotients, False)
	right_sequence.append(val)
	while max_s - min_s > EPS:
		s = (min_s + max_s) / 2.0
		(success, val) = extend_arrays(N, n, s, X1, Sigma1, C1, err1, A, quotients, False)
		
		# Success if resulting X[N, n:N] has modulus < R:
		if success:
			max_s = s
		else:
			min_s = s

	s = (max_s + min_s) / 2.0
	extend_arrays(N, n, s, X1, Sigma1, C1, err1, A, quotients, True)

	print("ANS:", np.dot(X1[N], C1[N]))

	return (X1, err1)



# =========================================================================================

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


# ============================================================================

def initial_state(N, M, X0):
	global P, R2

	X = np.zeros((N + 1, N + 1), dtype='double')
	for j in range(M):
		for i in range(M):
			X[i][j] = X0[i][j]
		for i in range(M, N + 1):
			X[i][j] = X0[N][j]

	for j in range(M, N):
		for i in range(j + 1, N + 1):
			X[i][j] = np.sqrt(R2 / (N - M))
	return X


def initial_guess(N):
	X = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(i):
			# X[i][j] = 1.0/np.sqrt(j + 1)
			# X[i][j] = 0.01
			# X[i][j] = 1.0/np.sqrt(j + 1) + random.random() / 10.0
			# X[i][j] = 0.071
			X[i][j] = 1.0 - random.random()
			# X[i][j] = (1.0 - 2.0 * random.random())/5.0
			# if j == 0:
			# 	X[i][j] = (-1)**i
			# elif j == 1:
			# 	X[i][j] = 0.0

	return X

def minimize():
	global FIXED_COORDS, modifiable, modifiableS, Z, R2
	FIXED_COORDS = False
	X0 = None
	S0 = None
	Z = None
	N = None
	M = None

	# Reading from input file
	# Indices 0 to M - 1 are meant for initial data.
	# Indices M to N - 1 are meant for evaluation points.
	# Index N is meant for output.
	if False:
		INPUT = read_from_file("in.txt")
		X0 = np.array(INPUT[0])
		N = len(X0) - 1
		Z0 = np.array(INPUT[2])
		M = len(Z0)

		modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
		# for i in range(M, N + 1):
		for i in range(M, N):
			for j in range(i):
				modifiable[i][j] = True

		modifiableS = [True for i in range(N)]
		for i in range(M):
			modifiableS[i] = False
		for i in range(M, N):
			modifiableS[i] = False
		S0 = get_S(N, INPUT[1])

		Z = np.eye(N + 1, dtype='double')
		for i in range(M):
			for j in range(M):
				Z[i][j] = Z0[i][j]
	else:
		N = 2
		M = 0
		P = 3.0

		Z = np.eye(N + 1, dtype='double')

		X0 = initial_guess(N)
		# X0[0][0] = 0.5798527
		# X0[0][1] = 0.32425464266
		# X0[0][2] = 0.163347653754
		# X0[1][0] = -0.57285235
		# X0[1][1] = -0.067735435654
		# X0[1][2] = -0.35463773
		# X0[2][0] = 0.867143785681
		# X0[2][1] = 0.538572989
		# X0[2][2] = -0.14325635
		# X0[N][0] = 0.0
		# X0[N][1] = 0.0
		# X0[N][2] = 0.0

		# X0[0][0] = -0.3
		# X0[0][1] = 0.0
		# X0[1][0] = 0.3
		# X0[1][1] = 0.0
		# X0[N][0] = 0.0
		# X0[N][1] = 0.0

		S0 = np.ones((N), dtype='double')
		# S0[0] = 4.418927941
		# S0[1] = 0.834729724
		# S0[2] = 2.943724987

		# S0[0] = 0.0
		# S0[1] = 0.0

		modifiableS = [True for i in range(N)]
		for i in range(M):
			modifiableS[i] = False
		for i in range(M, N):
			modifiableS[i] = False
			# S0[i] = (N - M)/P
			# S0[i] = (N - M)/P * 10**i
			S0[i] = np.random.random()
			S0[i] = S0[i] / (1.0 - S0[i])

		# Z[1][0] = -0.58743289793
		# Z[2][0] = 0.252353
		# Z[1][2] = 1.343151533
		# Z[2][1] = -0.375734553
		# Z[0][1] = 0.83925782935
		# Z[0][2] = -0.2342352359

		# Z[0][0] = 0.0
		# Z[0][1] = -2.0


		modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]
		for i in range(M, N):
		# for i in range(M, N + 1):
			for j in range(i):
				modifiable[i][j] = True
		for j in range(M):
			modifiable[N][j] = False
		# for i in range(N + 1):
		# 	modifiable[i][1] = False
		# 	X0[i][1] = 0.0

		# for i in range(3, N, 2):
		# 	modifiable[i][i - 1] = False
		# 	X0[i][i - 1] = 0.0
		
		# R2 = 6.8967142387450595
		# X0 = initial_state(N, M, X0)

	# print(X0)
	# print()
	# err0 = get_err2(N, S0)
	# while True:
	# # for i in range(2):
	# 	(X0, err0) = iterate(N, M, X0, err0)
	# 	time.sleep(0.1)
	# 	# print()
	# return
		
	
	fun = lambda v: calc(get_args(v, N, M, X0, S0))

	# (ANS0, grad) = fun(get_v(N, X0, S0))
	# for i in range(N - 1):
	# 	X1 = initial_guess(N)
	# 	X2 = initial_guess(N)
	# 	S1 = np.copy(S0)
	# 	S2 = np.copy(S0)
	# 	# X1[N - 1][i] += 0.001
	# 	# X2[N - 1][i] -= 0.001
	# 	X1[N][i] += 0.001
	# 	X2[N][i] -= 0.001
	# 	# S1[M + i] += 0.001
	# 	# S2[M + i] -= 0.001

	# 	(ANS1, _1) = fun(get_v(N, X1, S1))
	# 	(ANS2, _2) = fun(get_v(N, X2, S2))
	# 	# print((ANS1 - ANS2) / 0.002, grad[((N - 2) * (N - 1) // 2) - (M * (M - 1) // 2) + i])
	# 	print((ANS1 - ANS2) / 0.002, grad[(N * (N - 1) // 2) - (M * (M - 1) // 2) + i])
	# 	# print((ANS1 - ANS2) / 0.002, grad[(N * (N + 1) // 2) - (M * (M - 1) // 2) + i])


	# fun(get_v(N, X0, S0))
	# for j in range(M):
	# 	modifiable[N][j] = False
	RES = scipy.optimize.minimize(fun, get_v(N, X0, S0), jac=True, options={"gtol": 0.0000000000000000000001}, method='BFGS')
	
	print("=================================================================")
	print(RES)

	(N, M, X1, S1) = get_args(RES.x, N, M, X0, S0)
	print(X1)
	print(get_err2(N, S1))

	with open("out.txt", "a") as f:
		f.write(str(X1))
		f.write("\n")
		f.write(str(get_err2(N, S1)))
		f.write("\n")
		f.write(str(Z[:M, :M]))
		f.write("\n\n\n")

	# return (fun(RES.x)[0], X1)


minimize()

# best_ans = 100000.0
# best_X = None

# for k in range(10):
# 	print(k)
# 	(ans, X) = minimize()
# 	print(ans)

# 	if ans < best_ans:
# 		best_ans = ans
# 		best_X = X

# print(best_ans)
# print(best_X)
