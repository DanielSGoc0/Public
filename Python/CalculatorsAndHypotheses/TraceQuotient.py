# Check if one can bound A D A / tr(D A D A)
# for D positive diagonal and A positive definite.
# Actually, a bunch of other stuff as well...
import numpy as np
from scipy.stats import ortho_group

# N = 2

# BEST = 0.0
# BEST_K = None
# BEST_D1 = None
# BEST_D2 = None
# BEST_K1 = None
# BEST_K2 = None

# for attempt in range(100000):
# 	K = (np.diag(np.random.random((N))) + 99.0) / 100.0
# 	U = ortho_group.rvs(dim=N)
# 	K = U @ K @ U.T
	
# 	D1 = np.diag(np.random.random((N)))
# 	D2 = np.diag(np.random.random((N)))
# 	for i in range(N):
# 		D1[i][i] = (D1[i][i] / (1.0 - D1[i][i]))**3
# 		# D2[i][i] = D2[i][i] / (1.0 - D2[i][i])

# 		D2[i][i] = D1[i][i] * np.exp(np.random.normal(0.0, 0.01, (1)))
# 	K1 = K - K @ np.linalg.solve(K + D1, K)
# 	K2 = K - K @ np.linalg.solve(K + D2, K)

# 	linf = 0.0
# 	l2 = 0.0
# 	for i in range(N):
# 		linf = max(linf, np.abs(K1[i][i] - K2[i][i]))
# 		l2 += (K1[i][i] - K2[i][i])**2

# 	RES = np.abs(K1[1][0] - K2[1][0]) / linf

# 	if RES > BEST:
# 		BEST = RES
# 		BEST_K = K
# 		BEST_D1 = D1
# 		BEST_D2 = D2
# 		BEST_K1 = K1
# 		BEST_K2 = K2
# 	print(BEST, '\t', RES)

# print()
# print(BEST)
# print(BEST_K)
# print(BEST_D1)
# print(BEST_D2)
# print(BEST_K1)
# print(BEST_K2)

def distance(N, K, D1, D2):
	D1inv = np.zeros((N, N), dtype='double')
	D2inv = np.zeros((N, N), dtype='double')
	for i in range(N):
		D1inv[i][i] = 1.0 / D1[i][i]
		D2inv[i][i] = 1.0 / D2[i][i]

	K1 = K - K @ np.linalg.solve(K + D1inv, K)
	K2 = K - K @ np.linalg.solve(K + D2inv, K)

	return np.sqrt(-np.trace((D1 - D2) @ (K1 - K2)))

def conditioned(N, K, D):
	Dinv = np.zeros((N, N), dtype='double')
	for i in range(N):
		Dinv[i][i] = 1.0 / D[i][i]

	return K - K @ np.linalg.solve(K + Dinv, K)

def trln(N, K, D):
	det = np.linalg.det(np.eye(N, dtype='double') + D @ K)
	return np.log(det)


SUP = 0.0
BEST_D1 = None
BEST_D2 = None
BEST_D3 = None
BEST_K = None
BEST_K1 = None
BEST_K2 = None

N = 2
for attempt in range(10000):
	K = np.diag(np.random.random((N)))
	p = np.random.random()
	p = p / (1.0 - p)
	# p = 0.0
	for i in range(N):
		K[i][i] = (p + K[i][i])/(1.0 + p)
	# for i in range(N):
	# 	K[i][i] = K[i][i] / (1.0 - K[i][i])
	U = ortho_group.rvs(dim=N)
	K = U @ K @ U.T

	# K = np.eye(N, dtype='double')

	p = np.random.random()
	p = p / (1.0 - p)
	# p = 0.0
	D = np.diag(np.random.random((N)))
	D1 = np.diag(np.random.random((N)))
	D2 = np.diag(np.random.random((N)))
	D3 = np.diag(np.random.random((N)))
	for i in range(N):
		D[i][i] = (D[i][i] / (1.0 - D[i][i]))**2
		D1[i][i] = (D1[i][i] / (1.0 - D1[i][i]))**2
		D1[i][i] = (p*D[i][i] + D1[i][i])/(p + 1.0)
		D2[i][i] = (D2[i][i] / (1.0 - D2[i][i]))**2
		D2[i][i] = (p*D[i][i] + D2[i][i])/(p + 1.0)
		D3[i][i] = (D3[i][i] / (1.0 - D3[i][i]))**2
		D3[i][i] = (p*D[i][i] + D3[i][i])/(p + 1.0)
		# D2[i][i] = D1[i][i] * np.exp(np.random.normal(0.0, 0.01, (1)))
		# D3[i][i] = D1[i][i] * np.exp(np.random.normal(0.0, 0.01, (1)))

	# K1 = conditioned(N, K, D1)
	# K2 = conditioned(N, K, D2)
	# K1D = conditioned(N, K, D1 + D)
	# K2D = conditioned(N, K, D2 + D)

	# print(K1 - K2 - (K1 @ (D2 - D1) @ K2))

	# dist = -np.trace((D1 - D2) @ (K1 - K2))
	# dist2 = np.dot(np.diag(D1 - D2), np.multiply(K1, K2) @ np.diag(D1 - D2))
	# print(dist - dist2)
	# distD = -np.trace((D1 - D2) @ (K1D - K2D))
	# print(dist / distD)
	# dist = 0.0
	# distD = 0.0
	# for i in range(N):
	# 	for j in range(N):
	# 		# dist += (K1[i][j] - K2[i][j]) * (D1[i][i] - D2[j][j])
	# 		# distD += (K1D[i][j] - K2D[j][i]) * (D1[i][i] - D2[j][j])
	# 		# dist += K1[0][i] * K[i][j] * K2[j][0] * (D1[i][i] - D2[j][j])
	# 		# distD += K1D[0][i] * K[i][j] * K2D[j][0] * (D1[i][i] - D2[j][j])
	# 		# dist += (K1[i][j] * K2[j][i]) * D1[i][i] * D2[j][j]
	# 		# distD += (K1D[i][j] * K2D[j][i]) * (D1[i][i] + D[i][i]) * (D2[j][j] + D[j][j])
	# 		dist += (K1[i][j] * K2[j][i]) * (D1[i][i] - D2[i][i]) * D2[j][j]
	# 		distD += (K1D[i][j] * K2D[j][i]) * (D1[i][i] + D[i][i]) * (D2[j][j] + D[j][j])
	# print(dist, distD)
	# # dist = 0.0
	# # for i in range(N):
	# # 	dist -= (D1[i][i] - D2[i][i]) * (K1[i][i] - K2[i][i])
	# # 	# dist -= (D1[i][i] - D2[i][i]) * np.log(K1[i][i] / K2[i][i])
	# RES = 0.0
	# for i in range(N):
	# 	for j in range(N):
	# 		# if i == j:
	# 		# 	continue
	# 		# RES = max(RES, (K1[i][j] - K2[i][j])**2 / (K[i][i] * K[j][j] * dist))
	# 		# RES = max(RES, abs(np.log(K1[i][j] / K2[i][j])) / dist)
	# 		RES = max(RES, (K1[i][j] - K2[i][j])**2 / (K1[i][i] * K2[j][j] * dist))
	# 		# RES = max(RES, (K1[i][j] - K2[i][j])**2 / (max(K1[i][i], K1[j][j]) * max(K2[i][i], K2[j][j]) * dist))
	# 		# RES = max(RES, (K1[i][j] - K2[i][j])**2 / (min(K1[i][i], K1[j][j]) * min(K2[i][i], K2[j][j]) * dist))
	# 		# RES = max(RES, (K1[i][j] - K2[i][j])**2 / (np.sqrt(K1[i][i] * K1[j][j] *  K2[i][i] * K2[j][j]) * dist))
	# 		# RES = max(RES, (K1[i][j] - K2[i][j])**2 / ((K1[i][i] * K1[j][j] / 2 + K2[i][i] * K2[j][j] / 2) * dist))
	# 		# RES = max(RES, (K1[i][j] - K2[i][j])**2 / (np.sqrt((K1[i][i] * K1[j][j])**2 / 2 + (K2[i][i] * K2[j][j])**2 / 2) * dist))

	# 		# RES = max(RES, (K1[i][j] - K2[i][j])**2 / (((K1[i][i] - K2[i][i])**2 + (K1[j][j] - K2[j][j])**2)))

	# # print(dist)
	
	# # RES = np.max(np.abs(K1 - K2)) / np.sqrt(dist)
	# # RES = np.max(np.abs(K1 - K2)) / np.max(np.abs(np.diag(K1 - K2)))
	# # RES = np.abs(K1[0][1] - K2[0][1]) / np.sqrt(dist)

	d12 = distance(N, K, D1, D2)
	d23 = distance(N, K, D2, D3)
	d13 = distance(N, K, D1, D3)

	# d12 = d12 / (1.0 + d12)
	# d13 = d13 / (1.0 + d13)
	# d23 = d23 / (1.0 + d23)

	# print(d12, d13, d23)
	# RES = d12 / (d13 + d23)
	# if RES > 0.999:
	# 	print(d12, d13, d23)
	# 	break

	# print(trln(N, K, D1))
	# print(trln(N, K, D2))
	# print(trln(N, K, D3))

	RES = (d12 + 1.0) / ((d13 + 1.0) * (d23 + 1.0))

	if SUP < RES:
		SUP = RES
		BEST_D1 = D1
		BEST_D2 = D2
		BEST_D3 = D3
		BEST_K = K
		# BEST_K1 = K1
		# BEST_K2 = K2
	print(SUP, '\t', RES)

print()
print(SUP)
print(BEST_K)
print(BEST_D1)
print(BEST_D2)
print(BEST_D3)
print(BEST_K1)
print(BEST_K2)
