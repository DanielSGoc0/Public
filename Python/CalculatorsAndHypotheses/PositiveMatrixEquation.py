# ???
# This was helpful for something, but I can't remember what
import numpy as np

X = np.random.random((3, 3))
# M = X @ X.T
M = np.zeros((3, 3), dtype='double')
for i in range(3):
	for j in range(3):
		M[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)

print(M)
P = 1.0

a = M[2][0] - M[2][1]
b = M[0][0] * M[2][1] - M[0][1] * M[2][0]
c = M[1][0] * M[2][1] - M[1][1] * M[2][0] - M[0][0] * M[2][1] + M[0][1] * M[2][0]

d1 = c / (a - P * b)
d0 = c / (P * (c + b) - a)

print(d0)
print(d1)

# print((a - P * b) / c)
# print(P - (a - P * b) / c)
# print(M[2][0] - M[2][1])
# print(M[0][0] * M[2][1] - M[0][1] * M[2][0])
# print(M[1][0] * M[2][1] - M[1][1] * M[2][0] - M[0][0] * M[2][1] + M[0][1] * M[2][0])

# d0 = np.random.random()
# d1 = np.random.random()

INV = np.linalg.inv(M + np.diag([d0, d1, 0.0]))
print(d0 * INV[2][0], d1 * INV[2][1])
