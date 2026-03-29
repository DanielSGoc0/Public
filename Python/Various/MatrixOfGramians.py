import numpy as np

N = 5
M = 20

X = np.random.normal(loc=0.0, scale=1.0, size=(N, M))
for i in range(N):
	s = np.dot(X[i], X[i])
	X[i] /= np.sqrt(s)

G = np.zeros((N, N), dtype='double')
for i in range(N):
	for j in range(N):
		G[i][j] = np.dot(X[i], X[j])

M = np.zeros((N, N), dtype='double')
for i in range(N):
	for j in range(N):
		a = min(i, j)
		b = max(i, j) + 1
		M[i][j] = np.linalg.det(G[a:b, a:b])
	
print(M)
print(np.linalg.det(M))
