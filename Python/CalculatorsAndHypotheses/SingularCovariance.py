# Let K be a singular covariance matrix.
# Is it true that the optimal matrix D need
# not be supported on all indices?
import numpy as np
from scipy.stats import ortho_group
import cvxpy as cp

P = 100.0
N = 2
for attempt in range(10000):
	print(attempt)

	K = np.diag(np.random.random((N)))
	for i in range(N):
		K[i][i] = K[i][i] * 0.99 + 0.5 * 0.01
		K[i][i] = K[i][i]/(1.0 - K[i][i])
	K[0][0] = 0.0
	U = ortho_group.rvs(dim=N)
	K = U @ K @ U.T

	KN = np.random.normal(0.0, 1.0, (N))
	# We must ensure that KN is perpendicular to U[:, 0]!
	KN -= (KN @ U[:, 0]) * U[:, 0]

	v1 = cp.Variable(N)
	expr1 = 2 * KN @ v1 - cp.quad_form(v1, K, assume_PSD=True) - (1.0/P) * cp.norm1(v1)**2
	objective1 = cp.Maximize(expr1)
	# v2 = cp.Variable(N)
	# objective2 = cp.Maximize(
	# 	2 * KN @ v2 - cp.quad_form(v2, K) - (0.5/P) * cp.norm1(v2)**2
	# )

	problem1 = cp.Problem(objective1)
	problem1.solve()
	# problem2 = cp.Problem(objective2)
	# problem2.solve()

	# print(problem.value)
	# print(v.value)

	AMT_NON_ZERO_1 = 0
	for val in v1.value:
		if abs(val) > 0.000000001:
			AMT_NON_ZERO_1 += 1

	# AMT_NON_ZERO_2 = 0
	# for val in v2.value:
	# 	if abs(val) > 0.000000001:
	# 		AMT_NON_ZERO_2 += 1

	# print(AMT_NON_ZERO_1, AMT_NON_ZERO_2)

	if AMT_NON_ZERO_1 == N:
	# if AMT_NON_ZERO_1 == N and AMT_NON_ZERO_2 < N:
		print("ok")

		print(v1.value)
		new_v1 = v1.value + 0.1 * U[:, 0]
		print(2 * KN @ new_v1 - new_v1 @ K @ new_v1 - (1.0/P) * np.linalg.norm(new_v1, 1)**2)
		new_v1 = v1.value
		print(2 * KN @ new_v1 - new_v1 @ K @ new_v1 - (1.0/P) * np.linalg.norm(new_v1, 1)**2)
		new_v1 = v1.value - 0.1 * U[:, 0]
		print(2 * KN @ new_v1 - new_v1 @ K @ new_v1 - (1.0/P) * np.linalg.norm(new_v1, 1)**2)

		print(K)
		print(KN)
		# print(problem1.value)
		# print(v1.value)
		break
