# Check the following hypothesis:
# If phi is an RBF in infinite dimensions,
# and X0, ... , XN are arbitrary points,
# then the function Sigma[N, :N] Sigma[:N, :N]^{-1} Sigma[:N, N]
# is decrasing as we increase distance of XN from span(Xi).
import numpy as np


def phi(r2):
	return 1.0 + 5.0 * np.exp(-r2)
	# return 1.0 + np.exp(-r2/2.0)
	# return np.exp(-r2/4.0) + 4.0 * np.exp(-r2 * 4.0)


N = 4
# for k in range(1000000):
# 	print(k)
# 	X = np.random.random((N + 1, N + 1))

# 	Sigma = np.zeros((N + 1, N + 1), dtype='double')
# 	for i in range(N):
# 		for j in range(N):
# 			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))
# 		Sigma[i][i] += 1.0
# 	INV = np.linalg.inv(Sigma[:N, :N])

# 	# A = 1.0 - 2.0 * np.random.random((N, N))
# 	# INV = A @ A.T

# 	vals = []
# 	for s in range(0, 10):
# 		# R2 = 1.3**s
# 		R2 = s
# 		for i in range(N + 1):
# 			Sigma[N][i] = phi(np.dot(X[i] - X[N], X[i] - X[N]) + R2)
# 			Sigma[i][N] = Sigma[N][i]

# 		vals.append(Sigma[N, :N] @ INV @ Sigma[:N, N])
	
# 	S = len(vals)
# 	diffs = [[0.0 for j in range(S)] for i in range(S)]
# 	for j in range(S):
# 		diffs[0][j] = vals[j]
# 	for i in range(1, S):
# 		for j in range(S - i):
# 			diffs[i][j] = diffs[i - 1][j] - diffs[i - 1][j + 1]
# 			if diffs[i][j] < 0.0:
# 				print(X)
# 				print(diffs)
# 				exit()
# 	# print(diffs)

for k in range(1000000):
	X = 0.5 * np.random.random((N + 1, N + 1))

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N):
		for j in range(N):
			Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))
	INV = np.linalg.inv(Sigma[:N, :N])
	
	u = np.array([np.exp(-np.dot(X[i] - X[N], X[i] - X[N]) / 2.0) for i in range(N)])
	# u = Sigma[N, :N]
	# u = np.random.random((N))
	alpha = np.random.random()
	alpha = (1.0 - alpha) / alpha
	beta = np.random.random()
	beta = (1.0 - beta) / beta
	v = np.power(u, alpha)
	w = np.power(u, beta)

	ANS = v @ INV @ w
	if ANS < 0.0:
		print("oops!", ANS)
		print(u)
		print(alpha)
		print(beta)
		print(v)
		print(w)
		print(INV)
		break
