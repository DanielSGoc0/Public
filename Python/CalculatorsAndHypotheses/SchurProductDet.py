# Let K be positive definite.
# Let D be positive semi-definite, and diagonal?
# Then consider A = (K + D)^-1 K and also B = A o A
# Must det(B) > 0?
import numpy as np
from scipy.stats import ortho_group


# T = np.array([[0.0, -1.0], [1.0, 0.0]])
# print(np.linalg.eigvals(T))

def generateA(f):
	K = np.diag(np.random.random((N)))
	for i in range(N):
		K[i][i] = K[i][i] / (1.0 - K[i][i])
	U = ortho_group.rvs(dim=N)
	K = U @ K @ U.T
	D = np.diag(np.random.random((N)))
	for i in range(N):
		D[i][i] = f * D[i][i] / (1.0 - D[i][i])
	U = ortho_group.rvs(dim=N)
	D = U @ D @ U.T
	A = D @ np.linalg.solve(K + D, K)
	# A = D @ K

	return A


for attempt in range(100000):
	if attempt % 100 == 0:
		print(attempt)

	N = 2

	A1 = generateA(1.0)
	A2 = generateA(1.0)
	# A3 = generateA(1.0)
	# # A2 = U @ A2 @ U.T
	# M = np.diag(np.random.random((N)))
	# U = ortho_group.rvs(dim=N)
	# M = U @ M @ U.T


	# B = A1
	B = np.multiply(A1, A2)
	# B = (M + A1) / 2.0
	# B = np.multiply(M, A1)
	# B = np.multiply(A, A2)
	# B = np.multiply(B, A1)
	# B = np.multiply(B, A2)
	# B = np.multiply(B, A3)
	# B = np.reciprocal(A)

	det = np.linalg.det(B)
	evals = np.linalg.eigvals(B)
	# det = np.linalg.det(A)
	# evals = np.linalg.eigvals(A)

	FAIL = False
	# if det < 0.0:
	# 	FAIL = True

	for v in evals:
		# if not np.isreal(v) or v < -0.000001 or v > 1.0000001:
		if not np.isreal(v) or v < -0.000001:
			FAIL = True

	# if FAIL:
	if FAIL or attempt == 0:
		if FAIL:
			print("oh")
		# print(K)
		# print(D)
		# print(A1)
		# print(A2)
		print(B)
		print(evals)
		print(det)
		if FAIL:
			break
