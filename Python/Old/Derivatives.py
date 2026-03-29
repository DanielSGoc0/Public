# Calculates the limiting behavior of fixed coefficients methods, when coefficients are
# scaled by some factor h going to zero.

import scipy
import numpy as np

# we consider points 0 = x0, ... , x{N-1}
N = 5

h = 0.012

X = np.zeros((N + 1, N), dtype='double')
DF = np.zeros((N, N), dtype='double')
DF[0][0] = 1.0
Sigma = [[np.eye(N, dtype='double') for j in range(N)] for i in range(N)]
sigma = np.eye(N, dtype='double')

for k in range(1, N):
	# at the beginning we assume the knowledge of:
	# sigma[:k, :k], Sigma[:k, :k], X[:k], DF[:k].
	# The rest we have to update.

	# update DF[k - 1][k - 1] and X[k]
	X[k] = X[k - 1] + h*DF[k - 1]

	# update sigma at positions k.
	for l in range(k):
		sigma[l][k] = np.exp(-np.dot(X[l] - X[k], X[l] - X[k]) / 2.0)
		sigma[k][l] = sigma[l][k]

	# update Sigma at positions k.
	for l in range(k):
		for i in range(k + 1):
			for j in range(k + 1):
				Sigma[l][k][i][j] -= (X[l][i] - X[k][i]) * (X[l][j] - X[k][j])
		Sigma[l][k] *= sigma[l][k]
		Sigma[k][l] = Sigma[l][k]

	# Finally, calculate the projection of DF[k] onto (e0, ... , e{k-1}).
	A = np.zeros((k*k, k*k), dtype='double')
	B = np.zeros((k, k*k), dtype='double')
	v = np.zeros((k*k), dtype='double')

	for i in range(k):
		for j in range(k):
			for l in range(k):
				for m in range(k):
					A[k*l + i][k*m + j] = Sigma[l][m][i][j]

	for i in range(k):
		for j in range(k):
			for m in range(k):
				B[i][k*m + j] = Sigma[k][m][i][j]

	for i in range(k):
		for l in range(k):
			v[k*l + i] = DF[l][i]

	res = np.dot(B, np.linalg.solve(A, v))
	for i in range(k):
		DF[k][i] = res[i]

	c = np.dot(sigma[k, :k], np.linalg.solve(sigma[:k, :k], sigma[k, :k]))
	DF[k][k] = np.sqrt(1 - c)

X[N] = X[N - 1] + h*DF[N - 1]

# print(X)
# print(DF)
# print(sigma)
# for l in range(k):
#     for m in range(k):
#         print(Sigma[l][m])

# print("======================= s = " + str(s) + " ========================")
# print(X[1] / h)
# print((X[2] - 2*X[1]) / h**2)
# print((X[3] - 3*X[2] + 3*X[1]) / h**3)
# print((X[4] - 4*X[3] + 6*X[2] - 4*X[1]) / h**4)
# print((X[5] - 5*X[4] + 10*X[3] - 10*X[2] + 5*X[1]) / h**5)
# X[N] = X[N - 1] + h*DF[N - 1]
for k in range(N + 1):
    print("===============  k = " + str(k) + "  ====================")
    for i in range(N):
        print(X[k][i] / h**(i + 1))
