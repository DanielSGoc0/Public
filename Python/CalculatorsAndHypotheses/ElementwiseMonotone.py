# Checks if matrix inequalities are preserved
# under certain element-wise operations
import numpy as np

N = 3
EPS = 0.000001
M = np.ones((N, N), dtype='double')
for i in range(N):
	M[i][i] -= N


def cond_less_than(X, Y):
	global EPS, N, M

	evals = np.linalg.eigvalsh(M @ (Y - X) @ M)
	print(evals)

	RES = 1
	for v in evals:
		# if v < EPS:
		# 	RES = min(RES, 0)
		if v < -EPS:
			RES = min(RES, -1)

	return RES

def less_than(X, Y):
	global EPS, N, M

	evals = np.linalg.eigvalsh(Y - X)
	# print(evals)

	RES = 1
	for v in evals:
		if v < EPS:
			RES = min(RES, 0)
		if v < -EPS:
			RES = min(RES, -1)

	return RES
	
def phi(X):
	Y = np.copy(X)
	n = len(Y)
	
	for i in range(n):
		for j in range(n):
			Y[i][j] = np.exp(Y[i][j] / 2.0)
	
	return Y



for k in range(100000):

	A = np.random.random((N, N))
	DA = [[np.dot(A[i] - A[j], A[i] - A[j]) for j in range(N)] for i in range(N)]
	DA = np.array(DA)

	B = np.random.random((N, N))
	DB = [[np.dot(B[i] - B[j], B[i] - B[j]) for j in range(N)] for i in range(N)]
	DB = np.array(DB)

	RES1 = cond_less_than(-DA, -DB)
	RES2 = cond_less_than(phi(-DA), phi(-DB))
	if RES1 * RES2 == -1:
		print("BAD")
		print(DA)
		print(DB)
		break
	else:
		print(RES1, RES2)
