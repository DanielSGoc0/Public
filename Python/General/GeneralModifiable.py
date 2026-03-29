# Covariance is exp(-||x - y||^2 / 2)
# A most general implementation of Maximization problem/Fixed Point Problem.
import numpy as np
import scipy
import scipy.linalg
from scipy.stats import ortho_group

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=200000)

R2 = None
M = None
L = None
P = None
modifiable = None
PRINTING = False

def get_v(N, X, S, modifiable_X, modifiable_S, E):
	v = []

	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable_X[i][j]:
				v.append(X[i][j])

	for i in range(N):
		if modifiable_S[i]:
			v.append(S[i])

	return v


def get_X(N, v, X0, S0, modifiable_X, modifiable_S, E):
	X = np.zeros((N + 1, N + 1), dtype='double')

	k = 0
	for i in range(N + 1):
		for j in range(N + 1):
			if modifiable_X[i][j]:
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

	return (X, S, E)


# ==============================================================================================================================

# X is strictly lower diagonal
def calc(args):
	global R2, M, L, P, modifiable, PRINTING
	X = args[0]
	S = args[1]
	E = args[2]
	N = len(X) - 1
	X = np.array(X)

	err2 = np.zeros((N + 1), dtype='double')

	for i in range(L[0]):
		err2[i] = S[i]

	for k in range(len(L) - 2):
		for i in range(L[k], L[k + 1]):
			W = 0.0
			for i in range(L[k], L[k + 1]):
				W += np.exp(S[i])
			for i in range(L[k], L[k + 1]):
				err2[i] = 1.0/(np.exp(S[i]) * P[k] / W) 

	# for i in range(N):
	# 	err2[i] = S[i]


	X[N, M:] *= np.sqrt(R2 / np.dot(X[N, M:], X[N, M:]))

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
			# Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0) + 1.0
			# Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) * 4.0) + 4.0 * np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 4.0)
			# Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) * 8.0) + 1.0 * np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 8.0) + 1
		# M[i][i] += err2[i]
		Sigma[i][i] += err2[i]
	C = scipy.linalg.cholesky(Sigma + E, lower=True)
	# print(C)
	# print(scipy.linalg.cholesky(M[:N, :N], lower=False))


	ANS = 0.0
	for k in range(len(L) - 1):
		ANS += np.sqrt(np.dot(X[N, L[k]:L[k + 1]], X[N, L[k]:L[k + 1]]) * np.dot(C[N, L[k]:L[k + 1]], C[N, L[k]:L[k + 1]]))


	if PRINTING:
		print()
		print("ANS:", ANS)
		# h = np.sqrt(np.dot(X[N, OFFSET:], X[N, OFFSET:]) / (1.0 - np.dot(C[N, :OFFSET], C[N, :OFFSET])))
		# print("h:", h)

		f = np.zeros((len(L)), dtype='double')
		for k in range(len(L) - 1):
			f[k] = np.sqrt(np.dot(X[N][L[k]:L[k + 1]], X[N][L[k]:L[k + 1]]) / np.dot(C[N][L[k]:L[k + 1]], C[N][L[k]:L[k + 1]]))

		w = np.zeros((len(L) - 1), dtype='double')
		for k in range(len(L) - 2):
			w[k] = f[k] - f[k + 1]
		print("w:", np.array(w))

		# ANS2 = (w[1] + w[0]) * (Sigma[N][N] - Sigma[N, :OFFSET] @ np.linalg.solve(Sigma[:OFFSET, :OFFSET] + P[:OFFSET, :OFFSET], Sigma[:OFFSET, N]))
		# ANS2 -= w[1] * (Sigma[N][N] - Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N] + P[:N, :N], Sigma[:N, N]))
		# ANS2 -= w[0] * (Sigma[N][N] - Sigma[N, :L1] @ np.linalg.solve(Sigma[:L1, :L1] + P[:L1, :L1], Sigma[:L1, N]))
		# print("ANS2", ANS2)
		# print(X[1][0] / X[2][0])
		# print(Sigma[1][0] * Sigma[2][1] - (Sigma[1][1] - err2[1]) * Sigma[2][0])
		# print((Sigma[0][0] - err2[0]) * Sigma[2][1] - Sigma[0][1] * Sigma[2][0])
		# print("X:", np.abs(EVAL_X[OFFSET:N]))
		# print("C:", C[N, OFFSET:N])
		# print("p:", 1.0 / err2[OFFSET:N])
		# m = np.zeros((N - OFFSET), dtype='double')
		# for k in range(OFFSET, N):
		# 	Sigma_k = np.copy(Sigma[k, :])
		# 	Sigma_k[k] = 1.0
		# 	for i in range(k, N):
		# 		m[k - OFFSET] += w[i - OFFSET] * (Sigma[k][N] -  Sigma_k[:(i + 1)] @ np.linalg.inv(Sigma[:(i + 1), :(i + 1)]) @ Sigma[:(i + 1), N])**2
		# print("m:", m)
		# print("e:", err2)
		# print("p:", TOTAL1, TOTAL2)
		print(X)
		print(err2)
		# print(E)

		# for i in range(L1, N):
		# 	print(np.dot(X[i, OFFSET:L1], X[N, OFFSET:L1])**2 / (np.dot(X[i, OFFSET:L1], X[i, OFFSET:L1]) * np.dot(X[N, OFFSET:L1], X[N, OFFSET:L1])))

		b = np.linalg.solve(Sigma[:N, :N] + E[:N, :N], Sigma[:N, N])
		print(b)
		# V = np.zeros((N + 1, N + 1), dtype='double')
		# for i in range(N):
		# 	for j in range(N + 1):
		# 		if modifiable[i][j]:
		# 			# V[i][j] = -np.dot(b[:N], E[:N, j]) * X[i][j] - np.multiply(b[:N], Sigma[:N, j]) @ X[:N, j] + X[N][j] * Sigma[N][j]
		# 			V[i][j] = 2.0 * np.concatenate([b, [-1]]) @ (X[i][j] * Sigma[:, i] - np.multiply(Sigma[:, i], X[:, j])) * b[i]
		# 			# V[i][j] = b @ (X[i][j] * Sigma[:N, j] - np.multiply(Sigma[:N, j], X[:N, j]))
		# print(V)

		# print(Sigma[L1:N, L1:N] + np.diag(np.divide(E[L1:N, :N] @ b, b[L1:N])))
		# print(np.multiply(E[L1:N, :N] @ b, b[L1:N]))

	return ANS
	# return np.linalg.inv(Sigma[:N, :N] + E[:N, :N])
	# return Sigma + E



# =========================================================================================

def minimize(seed0, seed1):
	global R2, M, L, P, modifiable, PRINTING
	np.random.seed(seed0)

	R2 = 1.5
	M = 0
	L = [1, 3, 5]
	N = L[len(L) - 1]
	L.append(N + 1)
	P = [1.0, 1.0]

	X0 = np.zeros((N + 1, N + 1), dtype='double')
	modifiable = [[False for j in range(N + 1)] for i in range(N + 1)]

	for k in range(len(L) - 1):
		for i in range(L[k], L[k + 1]):
			for j in range(L[k]):
				if j >= M:
					modifiable[i][j] = True
				X0[i][j] = np.random.normal(0.0, 1.0, (1))
	for j in range(M):
		X0[N][j] = 0.0
	

	# S0 = np.ones((N), dtype='double')
	S0 = np.random.random((N))
	for i in range(L[0]):
		S0[i] = 1000000000.0

	modifiable_S = [False for i in range(N)]
	for i in range(M):
		modifiable_S[i] = False


	# E = np.diag(np.random.random((N + 1)))
	# # for i in range(N + 1):
	# # 	E[i][i] = E[i][i]/(1.0 - E[i][i])
	# U = ortho_group.rvs(dim=N+1, random_state=seed0)
	# E = U @ E @ U.T
	# for i in range(N + 1):
	# 	for j in range(N + 1):
	# 		if i < L1 or j < L1 or i == N or j == N:
	# 			E[i][j] = 0.0
	# 		# if i == N or j == N:
	# 		# 	E[i][j] = 0.0

	if seed0 != seed1:
		np.random.seed(seed1)
		for i in range(0, N + 1):
			for j in range(0, N + 1):
				if modifiable[i][j]:
					X0[i][j] = np.random.normal(0.0, 1.0, (1))


	E = np.zeros((N + 1, N + 1), dtype='double')
	fun = lambda v: -calc(get_X(N, v, X0, S0, modifiable, modifiable_S, E))

	# for i in range(L1, N):
	# 	for j in range(M, L1):
	# 		X0[i][j] -= 0.001
	# 		ANS0 = calc([X0, S0, Z, E], 2)
	# 		X0[i][j] += 0.002
	# 		ANS1 = calc([X0, S0, Z, E], 2)
	# 		X0[i][j] -= 0.001

	# 		print(i, j, "derivative:", (ANS1 - ANS0)/0.002)
	# 		ANS1 = calc([X0, S0, Z, E], 2)

	# return

	# RES = scipy.optimize.minimize(fun, get_v(N, X0, modifiable), options={"gtol": 0.00000001, 'ftol': 0.0000000000000000001}, method='SLSQP')
	RES = scipy.optimize.minimize(fun, get_v(N, X0, S0, modifiable, modifiable_S, E), options={"gtol": 0.0000000000000000000001})
	# RES = scipy.optimize.minimize(fun, get_v(N, X0, S0, modifiable, modifiable_S, Z))
	# RES = scipy.optimize.minimize(fun, get_v(N, X0, S0, modifiable, modifiable_S, Z), options={"gtol": 0.00000001})
	print("=================================================================")
	PRINTING = True
	fun(RES.x)

	with open("out.txt", "a") as f:
		f.write(str(get_X(N, RES.x, X0, S0, modifiable, modifiable_S, E)[0]))
		f.write("\n")
		f.write(str(get_X(N, RES.x, X0, S0, modifiable, modifiable_S, E)[1]))
		f.write("\n\n")

	return

	X1 = get_X(N, RES.x, X0, S0, modifiable, modifiable_S, Z, E)[0]
	CLOSE_EVALS = True
	ALL_PARALLEL = True
	POSITIVE = True
	for i in range(L1, N):
		if np.dot(X1[i] - X1[N], X1[i] - X1[N]) > 8.0:
			CLOSE_EVALS = False
		if np.dot(X1[i, M:L1], X1[N, M:L1])**2 / (np.dot(X1[i, M:L1], X1[i, M:L1]) * np.dot(X1[N, M:L1], X1[N, M:L1])) < 0.9:
			ALL_PARALLEL = False
	

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = np.exp(-np.dot(X1[i] - X1[j], X1[i] - X1[j]) / 2.0)
	
	b = np.linalg.solve(Sigma[:N, :N] + E[:N, :N], Sigma[:N, N])
	d = np.multiply(E[L1:N, :N] @ b, b[L1:N])
	for k in range(N - L1):
		if d[k] < -0.01:
			POSITIVE = False
			break
	
	if CLOSE_EVALS and (not ALL_PARALLEL): 
		# print("=================================================================")
		PRINTING = True
		fun(RES.x)
		PRINTING = False
	elif CLOSE_EVALS and not POSITIVE:
		PRINTING = True
		fun(RES.x)
		PRINTING = False
	elif ALL_PARALLEL and seed0 != seed1:
		print(-fun(RES.x))



	with open("out.txt", "a") as f:
		f.write(str(get_X(N, RES.x, X0, S0, modifiable, modifiable_S, Z, E)[0]))
		f.write("\n")
		f.write(str(get_X(N, RES.x, X0, S0, modifiable, modifiable_S, Z, E)[1]))
		f.write("\n\n")


minimize(0, 3)

# for seed in range(1000):
# 	print()
# 	print("SEED:", seed)
# 	minimize(seed, seed)

