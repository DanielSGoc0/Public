# Checks the hypothesis that conditioned
# covariance is concave in conditioning distribution
import numpy as np


def phi(x0, x1):
	return np.exp(-np.dot(x0 - x1, x0 - x1) / 2.0)

D = 2
n = 4
N = 7

# Indices 0:n are for first block
# Indices n:N are for second block
# Index N is for output
for k in range(10000):
	print(k)
	X = np.random.random((N + 1, D))

	Sigma = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(N + 1):
		for j in range(N + 1):
			Sigma[i][j] = phi(X[i], X[j])

	P0 = np.random.random(n)
	P2 = np.random.random(n)
	P1 = np.zeros((n), dtype='double')
	for i in range(n):
		P0[i] = P0[i] / (1.0 - P0[i])
		P2[i] = P2[i] / (1.0 - P2[i])
		P1[i] = (P0[i] + P2[i]) / 2.0
	
	D0 = np.zeros((N + 1, N + 1), dtype='double')
	D1 = np.zeros((N + 1, N + 1), dtype='double')
	D2 = np.zeros((N + 1, N + 1), dtype='double')
	for i in range(n):
		D0[i][i] = 1.0/P0[i]
		D1[i][i] = 1.0/P1[i]
		D2[i][i] = 1.0/P2[i]

	# VAL0 = Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N] + D0[:N, :N], Sigma[N, :N]) - Sigma[N, :n] @ np.linalg.solve(Sigma[:n, :n] + D0[:n, :n], Sigma[N, :n])
	# VAL1 = Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N] + D1[:N, :N], Sigma[N, :N]) - Sigma[N, :n] @ np.linalg.solve(Sigma[:n, :n] + D1[:n, :n], Sigma[N, :n])
	# VAL2 = Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N] + D2[:N, :N], Sigma[N, :N]) - Sigma[N, :n] @ np.linalg.solve(Sigma[:n, :n] + D2[:n, :n], Sigma[N, :n])
	# VAL0 = Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N] + D0[:N, :N], Sigma[N, :N])
	# VAL1 = Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N] + D1[:N, :N], Sigma[N, :N])
	# VAL2 = Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N] + D2[:N, :N], Sigma[N, :N])
	VAL0 = np.linalg.inv(Sigma[:N, :N] + D0[:N, :N])[0][1]
	VAL1 = np.linalg.inv(Sigma[:N, :N] + D1[:N, :N])[0][1]
	VAL2 = np.linalg.inv(Sigma[:N, :N] + D2[:N, :N])[0][1]

	# print(VAL0, VAL1, VAL2)
	# print(VAL1 - (VAL0 + VAL2) / 2.0)
	# if VAL1 - (VAL0 + VAL2) / 2.0 < 0.0:
	if VAL1 - (VAL0 + VAL2) / 2.0 > 0.0:
		print("ok")
		print(VAL0, VAL1, VAL2)
		print(VAL1 - (VAL0 + VAL2) / 2.0)
		break
	