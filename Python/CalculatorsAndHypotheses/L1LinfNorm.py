# Used to check if there exists an upper bound
# on the L1Linf norm on positive definite matrices
import numpy as np
from scipy.stats import ortho_group

# for k in range(1000000):
# 	N = 3
# 	A = 1.0 - 2.0 * np.random.random((N, N))
# 	A = A @ A.T

# 	# d = np.random.random()
# 	# D = np.eye(N) * d / (1.0 - d)
# 	D = np.diag(np.random.random((N)))
# 	for i in range(N):
# 		D[i][i] = D[i][i] / (1.0 - D[i][i]) 

# 	v = np.ones((N))
# 	for i in range(N):
# 		if np.random.random() > 0.5:
# 			v[i] = -v[i]

# 	c = np.linalg.inv(A + D) @ v
# 	LHS = np.sum(np.abs(c))
# 	# LHS = np.sum(np.abs(np.linalg.inv(A + D)))
# 	RHS = np.sum(1.0 / np.diag(D))



# 	SUM = 0
# 	for i in range(N):
# 		if v[i] * c[i] < 0.0:
# 			SUM += 1


# 	# print(k)
# 	# print(k, SUM)
# 	# print(k, LHS, RHS)
# 	# print(k, SUM, LHS, RHS)

# 	if LHS > RHS:
# 	# if SUM1 > SUM2:
# 		print("BAD")
# 		print(A)
# 		print(D)
# 		print(v)
# 		break

MAX = 0.0
for k in range(1000000):
	N = 2

	a = 1.0 - 2.0 * np.random.random((N))
	b = 1.0 - 2.0 * np.random.random((N))

	evals = np.random.random((N))
	o = ortho_group.rvs(dim=N)
	A = o @ np.diag(evals) @ o.T

	RHS = 0.0
	for i in range(N):
		RHS += abs(a[i] * b[i])
	LHS = a @ A @ b

	# if LHS > 2.0 * RHS:
	# 	print("damn")
	# 	print(evals)
	# 	print(A)
	# 	print(a)
	# 	print(b)
	# 	print(LHS, RHS)
	# 	break

	MAX = max(MAX, LHS / RHS)

	print(k, '\t', MAX)
