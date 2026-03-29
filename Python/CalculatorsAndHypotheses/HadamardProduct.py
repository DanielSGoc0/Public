# Check the following:
# If K is PSD and has eigenvalue 0,
# Then K Hadamard K is singular
import numpy as np
from scipy.stats import ortho_group


K = [[1, 1, 0, 1], [1, 3, 2, 1], [0, 2, 2, 0], [1, 1, 0, 0]]
K = np.array(K)
print(K)
N = len(K) - 1
# print(np.linalg.eigvalsh(K))

# K = np.multiply(K, K)
# print(np.linalg.eigvalsh(K))

e = 0.00001
D1 = np.diag(np.array([1, 0, 1, 0]))
D2 = np.diag(np.array([0, 1, 0, 0]))

I1 = []
I2 = []
for i in range(N + 1):	
	if D1[i][i] > e:
		I1.append(i)
	if D2[i][i] > e:
		I2.append(i)

A1 = K - K[np.ix_(range(N + 1), I1)] @ np.linalg.solve(K[np.ix_(I1, I1)] + np.linalg.inv(D1[np.ix_(I1, I1)]), K[np.ix_(I1, range(N + 1))])
A2 = K - K[np.ix_(range(N + 1), I2)] @ np.linalg.solve(K[np.ix_(I2, I2)] + np.linalg.inv(D2[np.ix_(I2, I2)]), K[np.ix_(I2, range(N + 1))])

print(A1)
print(A2)
