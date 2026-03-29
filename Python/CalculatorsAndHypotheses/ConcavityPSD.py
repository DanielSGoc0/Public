# Check if expression A1 (D1 - D2) A2 (D1 - D2) A3
# is positive semi-definite
import numpy as np
from scipy.stats import ortho_group


N = 5
for attempt in range(10000):
	Kinv = np.diag(np.random.random((N)))
	for i in range(N):
		Kinv[i][i] = Kinv[i][i]/(1.0 - Kinv[i][i])
	U = ortho_group.rvs(dim=N)
	Kinv = U @ Kinv @ U.T

	# p = np.random.random() * 2.0 + 1.0
	p = np.random.random()
	D1 = np.diag(np.random.random((N)))
	D2 = np.diag(np.random.random((N)))
	D3 = np.zeros((N, N), dtype='double')
	for i in range(N):
		D1[i][i] = D1[i][i]/(1.0 - D1[i][i])
		D2[i][i] = D2[i][i]/(1.0 - D2[i][i])
		# D3[i][i] = p * D1[i][i] + (1.0 - p) * D2[i][i]
		D3[i][i] = D1[i][i] + D2[i][i]

	K = np.linalg.inv(Kinv)
	A1 = np.linalg.inv(Kinv + D1)
	A2 = np.linalg.inv(Kinv + D2)
	A3 = np.linalg.inv(Kinv + D3)

	# print(A1 - A2 - A1 @ (D2 - D1) @ A2)
	# print(A1 - A3 - (1.0 - p) * A1 @ (D2 - D1) @ A3)
	# print(A2 - A3 + p * A2 @ (D2 - D1) @ A3)

	# print(A1 @ (D1 - D2) @ A2 @ (D1 - D2) @ A3 - A2 @ (D1 - D2) @ A3 @ (D1 - D2) @ A1)
	# M = A1 @ (D1 - D2) @ A2 @ (D1 - D2) @ A3
	# F = A3 - p * A1 - (1.0 - p) * A2
	# print(F + p * (1.0 - p) * M)
	# evals = np.linalg.eigvalsh((M + M.T) / 2.0)
	# evals = np.linalg.eigvalsh(-F)
	evals = np.diag(K + A3 - (A1 + A2))
	# print(evals)

	for v in evals:
		if v < -0.000001:
			print("=======================")
			print(np.diag(D3))
			print(evals)
			break

	# if A1[0][1] * A2[0][1] < -0.00001:
	# 	print("here")
	# 	print(A1)
	# 	print(A2)
	# 	break
