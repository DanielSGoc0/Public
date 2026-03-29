# Evaluates the maximal covariance between
# some off-grid point and all the possible points
# on the grid.
import numpy as np

a = -1.0
b = 1.0
S = 9
R2 = 1.0

# X = np.linspace(a, b, S)
Xp, Yp = np.meshgrid(np.linspace(a, b, S), np.linspace(a, b, S))
X = np.array([[Xp[i // S][i % S], Yp[i // S][i % S]] for i in range(S**2)])
# print(X)

N = len(X)

def phi(r2):
	# return np.exp(-r2 / 2.0)
	return np.exp(-r2) + np.exp(-r2 / 4.0)

Sigma = np.zeros((N + 1, N + 1), dtype='double')
for i in range(N):
	for j in range(N):
		Sigma[i][j] = phi(np.dot(X[i] - X[j], X[i] - X[j]))
	Sigma[i][N] = phi(np.dot(X[i], X[i]) + R2)
	Sigma[N][i] = Sigma[i][N]

print(phi(2.0 * R2))
print(np.linalg.solve(Sigma[:N, :N], Sigma[:N, N]))
print(Sigma[N, :N] @ np.linalg.solve(Sigma[:N, :N], Sigma[:N, N]))
	
