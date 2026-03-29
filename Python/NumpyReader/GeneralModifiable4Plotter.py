# Plotter2D specifically for visualizing optimal algorithms from GeneralModifiable4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from pathlib import Path
import scipy.linalg

np.set_printoptions(suppress=True, formatter={'float_kind':'{:5.8f}\t'.format}, linewidth=200000, threshold=np.inf)

def read_number(txt, cursor):
	x = 0
	negative = False
	after_decimal = -1
	while True:
		c = txt[cursor]
		if c == ' ' or c == '\t' or c == '\n' or c == ']':
			break
		elif c == '-':
			cursor += 1
			negative = True
		elif c == '.':
			cursor += 1
			after_decimal = 0
			continue
		else:
			cursor += 1
			x *= 10
			x += ord(c) - ord('0')
		if after_decimal >= 0:
			after_decimal += 1
	if after_decimal >= 0:
		x = float(x)
		x /= 10**after_decimal
	if negative:
		x = -x
	return (x, cursor)

def create_array(txt, cursor):
	RES = []
	while True:
		if cursor == len(txt):
			return (RES, cursor)
		c = txt[cursor]
		if c == ' ' or c == '\t' or c == '\n':
			cursor += 1
			continue
		elif c == '[':
			cursor += 1
			output = create_array(txt, cursor)
			RES.append(output[0])
			cursor = output[1]
		elif c == ']':
			cursor += 1
			return (RES, cursor)
		else:
			output = read_number(txt, cursor)
			RES.append(output[0])
			cursor = output[1]

def read_from_file(filename):
	txt = Path(filename).read_text()
	RES = create_array(txt, 0)[0]
	return RES


# =============================================================================================================================

# 1. full initial data, N = 100, n = 3, R2 = 0.5, P = 1
# 2. full initial data, N = 200, n = 3, R2 = 0.5, P = 1
# 3. full initial data, N = 200, n = 3, R2 = 0.5, P = 10
# 4. full initial data, naive method, N = 100, n = 3, R2 = 0.5, P = 1
# 5. full initial data, naive method, N = 100, n = 3, R2 = 0.5, P = 10
# 6. full initial data, N = 1000, n = 3, R2 = 0.1, P = 1000
# 7. full initial data, N = 1000, n = 3, R2 = 0.01, P = 1000
# 8. full initial data, N = 1000, n = 3, R2 = 0.001, P = 1000
# 9. full initial data, N = 1000, n = 3, R2 = 0.01, P = 100

INPUT = read_from_file("../ContinuousMethods/8.txt")

X = np.array(INPUT[0])
sigma2 = np.array(INPUT[1])
N = len(INPUT[0]) - 1
n = 3
R2 = np.dot(X[N, n:], X[N, n:])
print("R2 =", R2)

X_axis = [[None for j in range(N - n)] for i in range(N - n)]
Y_axis = [[None for j in range(N - n)] for i in range(N - n)]
Z_axis = [[None for j in range(N - n)] for i in range(N - n)]

TOTAL = 1.0 / sigma2[n]
for i in range(N - n):
	for j in range(N - n):
		X_axis[i][j] = TOTAL
		Y_axis[j][i] = TOTAL
	if i != N - n - 1:
		TOTAL += 1.0/sigma2[n + i + 1]


# ==============================================================================================================================

Sigma = np.zeros((N + 1, N + 1), dtype='double')
for i in range(N + 1):
	for j in range(N + 1):
		Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j]) / 2.0)
	Sigma[i][i] += sigma2[i]
C = scipy.linalg.cholesky(Sigma, lower=True)

Isigma = np.zeros((N + 1, N + 1), dtype='double')
Isigma2 = np.zeros((N + 1, N + 1), dtype='double')
Iexp = np.zeros((N + 1, N + 1), dtype='double')
for i in range(n, N):
	Isigma[i][i] = np.sqrt(sigma2[i])
	Isigma2[i][i] = sigma2[i]

Xsigma = X @ Isigma
Csigma = C @ Isigma


h = np.sqrt((1.0 - np.dot(C[N, :n], C[N, :n])) / R2)
print("h =", h)
W = np.zeros((N - n, N - n), dtype='double')
for i in range(n, N):
	for j in range(n, N):
		W[i - n][j - n] = h * C[i][j] - Sigma[N][i] * X[i][j]

w = np.zeros((N + 1), dtype='double')
w[n:N] = np.linalg.solve(W.T, X[N, n:N])
u = C.T @ w

K = np.zeros((N + 1, N + 1), dtype='double')
for i in range(N + 1):
	for j in range(N + 1):
		K[i][j] = C[N][max(i, j)] * u[min(i, j)]
M = scipy.linalg.solve_triangular(C.T, K.T, lower=False)
M = scipy.linalg.solve_triangular(C.T, M.T, lower=False)

A = np.multiply(M, Sigma)
# U satisfies U U* = -A and is upper triangular
# U = np.flip(scipy.linalg.cholesky(-np.flip(A[n:N, n:N], (0, 1)), lower=True), (0, 1))

wsigma2 = Isigma2 @ w
usigma = Isigma @ u
Msigma = Isigma2 @ M @ Isigma2

# =============================================================================================================================

# This is where we write all of the printing and experiments...

# print(np.diag(Msigma))
# print(usigma)

print(X[N][N]**2)
# print(np.divide(np.diag(X[(n + 1):N, n:(N - 1)]), X[N, n:(N - 1)])**2)
# print(np.divide(u[n:N], C[N, n:N]))
Z_axis = Xsigma[n:N, n:N]
# Z_axis = Msigma[n:N, n:N]
# Z_axis = Xsigma[n:N, n:N]
# Z_axis = Xsigma[(n + 1):, (n + 1):]
# for i in range(n, N):
# 	Sigma[i][i] -= sigma2[i]
# Z_axis = np.linalg.inv(Sigma[:N, :N])[n:N, n:N]
# Z_axis = np.log(Sigma[n:N, n:N])
# Z_axis = Msigma[n:N, n:N]
# print((h * usigma - Xsigma[N])[n:N])
# print(wsigma2)
# for i in range(n, N - 1):
# 	# print(C[N][i] / u[i])
# 	# print(C[N][i + 1] / u[i + 1] - C[N][i] / u[i])
# 	# print(u[i + 1] / C[N][i + 1] - u[i] / C[N][i])
# 	# print(Csigma[N][i + 1] * usigma[i] - Csigma[N][i] * usigma[i + 1])
# 	# print((Csigma[N][i + 1] * usigma[i] - Csigma[N][i] * usigma[i + 1])/usigma[i]**2)
# 	print((Csigma[N][i + 1] * usigma[i] - Csigma[N][i] * usigma[i + 1])/Csigma[N][i]**2)
# for j in range(n, N):
# 	for i in range(n, N):
# 		# if i - 1 > j:
# 		# 	Z_axis[i - n][j - n] = Isigma2[j][j] * (X[i][j]**2 - X[i - 1][j]**2)
# 		# if i - 2 > j:
# 		# 	Z_axis[i - n][j - n] = Isigma2[j][j] * (X[i][j]**2 - 2 * X[i - 1][j]**2 + X[i - 2][j]**2)
# 		if i - 3 > j:
# 			Z_axis[i - n][j - n] = Isigma2[j][j] * (X[i][j]**2 - 3 * X[i - 1][j]**2 + 3 * X[i - 2][j]**2 - X[i - 3][j]**2)
# 		else:
# 			Z_axis[i - n][j - n] = 0.0
# print(Z_axis)
# print(C[N] - h * X[N])


# =============================================================================================================================

for i in range(N - n):
	if i == 0:
		Z_axis[i][i] = Z_axis[i + 1][i]
	else:
		Z_axis[i][i] = Z_axis[i][i - 1]

# for i in range(N - n):
# 	Z_axis[i][0] = Z_axis[i][1]
# 	Z_axis[N - n - 1][i] = Z_axis[N - n - 2][i]
# Z_axis[0][0] = Z_axis[1][1]
# Z_axis[N - n - 1][N - n - 1] = Z_axis[N - n - 2][N - n - 2]

for i in range(0, N - n):
	for j in range(i, N - n):
		X_axis[i][j] = X_axis[j][j]
		Y_axis[i][j] = Y_axis[j][j]
		Z_axis[i][j] = Z_axis[j][j]

print(Z_axis[980:, 980:])

for i in range(N - n):
	for j in range(N - n):
		if j < 100 or j > N - n - 100:
			Z_axis[i][j] = 0.0


# mapping = [[(i, j) for j in range(N - n)] for i in range(N - n)]
# for i in range(N - n):
# 	for j in range(N - n):
# 		if j > i:
# 			mapping[i][j] = mapping[i][j - 1]
# 		else:
# 			x = min(N - n - 2, max(1, j))
# 			y = min(N - n - 2, max(1, i))
# 			while x > y - 2:
# 				if y < N - n - 2:
# 					y += 1
# 				else:
# 					x -= 1
# 			mapping[i][j] = (y, x)

# for i in range(N - n):
# 	for j in range(N - n):
# 		X_axis[i][j] = X_axis[mapping[i][j][0]][mapping[i][j][1]]
# 		Y_axis[i][j] = Y_axis[mapping[i][j][0]][mapping[i][j][1]]
# 		Z_axis[i][j] = Z_axis[mapping[i][j][0]][mapping[i][j][1]]


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.autoscale(enable=None, axis="x", tight=True)

ax.set_xlabel('used up precision')
ax.set_ylabel('coordinate index')
ax.set_zlabel('X')
surf = ax.plot_surface(X_axis, Y_axis, Z_axis, cmap=cm.coolwarm, linewidth=0, antialiased=False)

plt.show()
