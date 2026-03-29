# Plotter2D specifically for visualizing optimal algorithms
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

# ========================================================================

# INDEX = 0: no initial data, p = 1, N = 50, OFFSET = 0
# INDEX = 2: no initial data, p = 1, N = 100, OFFSET = 0
# INDEX = 4: no initial data, p = 10, N = 50, OFFSET = 0
# INDEX = 6: no initial data, p = 0.1, N = 50, OFFSET = 0
# INDEX = 8: full initial data with Z, p = 1, N = 100, OFFSET = 3
# INDEX = 10: full initial data with Z, p = 1, N = 50, OFFSET = 3
# INDEX = 12: full initial data with Z, p = 10, N = 100, OFFSET = 3
# INDEX = 14: full initial data with Z, p = 100, N = 100, OFFSET = 3
# INDEX = 16: + 2*e^(-|XN|^2 / 2), p = 1, N = 80, OFFSET = 0
# INDEX = 18: + 2*e^(-|XN|^2 / 2), p = 1, N = 40, OFFSET = 0
# INDEX = 20: + 2*e^(-|XN|^2 / 2), p = 2, N = 40, OFFSET = 0
# INDEX = 22: + 2*e^(-|XN|^2 / 2), p = 2, N = 80, OFFSET = 0
# INDEX = 24: full initial data with Z and fixed error, p = 1, N = 100, OFFSET = 3
# INDEX = 26: + 2*e^(-|XN|^2 / 2), p = 10, N = 100, OFFSET = 0
# INDEX = 28: + 2*e^(-|XN|^2 / 2) and fixed error, p = 1, N = 100, OFFSET = 0
# INDEX = 30: + 2*e^(-|XN|^2 / 2) and fixed error, p = 10, N = 100, OFFSET = 0

OFFSET = 0
INDEX = 0

INPUT = read_from_file("numpy_input.txt")
MAX = len(INPUT[INDEX]) - 1
MIN_X = 0.0
MAX_X = 1.0
MIN_Y = 0.0
MAX_Y = 1.0

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.autoscale(enable=None, axis="x", tight=True)

# Z = np.array([[f(X[i][j], Y[i][j]) for j in range(MAX + 1)] for i in range(MAX + 1)])
X = [[MIN_X + (MAX_X - MIN_X) * i / MAX for j in range(MAX - OFFSET)] for i in range(MAX - OFFSET)]
Y = [[MIN_Y + (MAX_Y - MIN_Y) * j / MAX for j in range(MAX - OFFSET)] for i in range(MAX - OFFSET)]
Z = np.array([[None for j in range(MAX - OFFSET)] for i in range(MAX - OFFSET)])
TOTAL = 1.0/INPUT[INDEX + 1][OFFSET]
for i in range(MAX - OFFSET):
	for j in range(MAX - OFFSET):
		X[i][j] = TOTAL
	if i != MAX - 1 - OFFSET:
		TOTAL += 1.0/INPUT[INDEX + 1][i + 1 + OFFSET]
for i in range(MAX - OFFSET):
	for j in range(MAX - OFFSET):
		Y[i][j] = X[j][i]

# ==============================================================================================================================

Xinput = np.array(INPUT[INDEX])
N = MAX
Sigma = np.zeros((N + 1, N + 1), dtype='double')
for i in range(N + 1):
	for j in range(N + 1):
		Sigma[i][j] = np.exp(-np.dot(Xinput[i] - Xinput[j], Xinput[i] - Xinput[j]) / 2.0)
	Sigma[i][i] += INPUT[INDEX + 1][i]
C = scipy.linalg.cholesky(Sigma, lower=True)
# for i in range(OFFSET, N + 1):
# 	Sigma[i][i] -= INPUT[INDEX + 1][i]

Isigma = np.zeros((N + 1, N + 1), dtype='double')
Isigma2 = np.zeros((N + 1, N + 1), dtype='double')
Iexp = np.zeros((N + 1, N + 1), dtype='double')
for i in range(OFFSET, N):
	Isigma[i][i] = np.sqrt(INPUT[INDEX + 1][i])
	Isigma2[i][i] = INPUT[INDEX + 1][i]
	Iexp[i][i] = np.exp(-np.dot(Xinput[i][OFFSET:], Xinput[i][OFFSET:]) / 2.0)
for i in range(OFFSET):
	Isigma[i][i] = 1.0


Xsigma = np.dot(Xinput, Isigma)
Csigma = np.dot(C, Isigma)

K = np.zeros((N + 1, N + 1), dtype='double')
Kprime = np.zeros((N + 1, N + 1), dtype='double')
for i in range(N + 1):
	for j in range(N + 1):
		K[i][j] = C[N][max(i, j)] * Xinput[N][min(i, j)]
		Kprime[i][j] = Xinput[N][max(i, j)] * C[N][min(i, j)]

M = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(K), lower=False)
M = scipy.linalg.solve_triangular(np.transpose(C), np.transpose(M), lower=False)
MSigma = np.zeros((N + 1, N + 1), dtype='double')
# Mchanged = np.zeros((N + 1, N + 1), dtype='double')
# Mchanged = Isigma @ Isigma @ M @ Isigma @ Isigma
m = M[N - 1][N - 1] * INPUT[INDEX + 1][N - 1]**2

for i in range(N + 1):
	for j in range(N + 1):
		MSigma[i][j] = M[i][j] * Sigma[i][j]
		# Mchanged[i][j] = Mchanged[i][j] / Sigma[i][j]**2
		# if i == j or i < OFFSET or j < OFFSET or i == N or j == N:
		# 	Mchanged[i][j] = 1.0
		# else:
		# 	Mchanged[i][j] = np.log(m/Mchanged[i][j])/np.log(Sigma[i][j])
		# 	# Mchanged[i][j] = -(Mchanged[i][j]/m - 1.0)/np.log(Sigma[i][j])

print("m:", m)
print(np.diag(M) * np.diag(Isigma2) * np.diag(Isigma2))

# Z = np.zeros((N - OFFSET, N - OFFSET), dtype='double')
# for i in range(N - OFFSET):
# 	for j in range(N - OFFSET):
# 		# if j <= i:
# 		# 	Z[i][j] = C[i + OFFSET + 1][j + OFFSET] / Xinput[i + OFFSET + 1][j + OFFSET]
# 		# Z[i][j] = Xinput[i + OFFSET][j + OFFSET] / Xinput[N][j + OFFSET]

# M = np.transpose(M)
# Z = M[(OFFSET + 1): , OFFSET:N]
# Z = (Isigma @ Isigma @ M @ Isigma @ Isigma)[OFFSET:N , OFFSET:N]
# Z = (Xinput @ Isigma)[(OFFSET + 1): , OFFSET:N]
# Z = (C @ Isigma)[(OFFSET + 1): , OFFSET:N]
# for i in range(N - OFFSET - 11):
# 	for j in range(N - OFFSET - 11):
# 		Z[i][j] = np.linalg.det(Xsigma[i:(i + 11):5, j:(j + 11):5])
# Z = Xinput[(OFFSET + 1): , OFFSET:N]
# Z = np.divide((1.0 - (Xinput @ Isigma))[(OFFSET + 1): , OFFSET:N], X[(OFFSET + 1): , OFFSET:N])
# Z = (C @ Isigma)[(OFFSET + 1): , OFFSET:N]
# Z = (Isigma @ np.linalg.inv(C) @ Isigma2)[OFFSET:N , OFFSET:N]
# for i in range(N - OFFSET):
	# for j in range(N - OFFSET):
		# Z[i][j] /= Xinput[N][j + OFFSET]
# for j in range(OFFSET, N):
# 	for i in range(OFFSET, N):
# 		# if i - 1 > j:
# 		# 	Z[i - OFFSET][j - OFFSET] = Isigma[j][j] * (Xinput[i][j] - Xinput[i - 1][j])
# 		# if i - 2 > j:
# 		# 	Z[i - OFFSET][j - OFFSET] = Isigma[j][j] * (Xinput[i][j] - 2 * Xinput[i - 1][j] + Xinput[i - 2][j])
# 		if i - 3 > j:
# 			Z[i - OFFSET][j - OFFSET] = Isigma[j][j] * (Xinput[i][j]**2 - 3 * Xinput[i - 1][j]**2 + 3 * Xinput[i - 2][j]**2 - Xinput[i - 3][j]**2)
# 		else:
# 			Z[i - OFFSET][j - OFFSET] = 0.0

# print(Xsigma[OFFSET:10, :10])
# print(100000*Z)

# calculates matrix inverse of a_{max(i, j)} * b_{min(i, j)}
def matrix_inverse(a, b, N):
	M = np.zeros((N, N), dtype='double')
	M[0][0] = b[1] / b[0] / (a[0] * b[1] - a[1] * b[0])
	M[N - 1][N - 1] = a[N - 2] / a[N - 1] / (a[N - 2] * b[N - 1] - a[N - 1] * b[N - 2])

	for k in range(0, N - 1):
		M[k][k + 1] = 1.0 / (a[k + 1] * b[k] - a[k] * b[k + 1])
		M[k + 1][k] = M[k][k + 1]
	
	for k in range(1, N - 1):
		M[k][k] = (a[k - 1] * b[k + 1] - a[k + 1] * b[k - 1]) / ((a[k] * b[k + 1] - a[k + 1] * b[k]) * (a[k - 1] * b[k] - a[k] * b[k - 1]))

	return M

# factor = C[N][N - 1] / Xinput[N][N - 1]
# Inner = matrix_inverse(Xinput[N, OFFSET:N], factor * Xinput[N, OFFSET:N] - C[N, OFFSET:N], N - OFFSET)
# Z = Xinput[(OFFSET+1):, OFFSET:N] @ Inner @ Xinput[(OFFSET+1):, OFFSET:N].T
# Z = Xinput[OFFSET:N, OFFSET:N] @ Inner @ Xinput[OFFSET:N, OFFSET:N].T
# print(np.diag(Z) / np.diag(Isigma2[(OFFSET + 1):, (OFFSET + 1):]))
# print(np.diag(Z))

# Inner = matrix_inverse(C[N, (OFFSET+1):], Xinput[N, (OFFSET+1):], N - OFFSET)
# Z = C[(OFFSET+1):, (OFFSET+1):] @ Inner @ C[(OFFSET+1):, (OFFSET+1):].T
# Inner = matrix_inverse(C[N, (OFFSET+1):], Xinput[N, (OFFSET+1):], N - OFFSET)
# Z = C[(OFFSET+1):, (OFFSET+1):] @ Inner @ C[(OFFSET+1):, (OFFSET+1):].T

# Z = np.dot(Isigma, np.dot(np.dot(Xinput[N], C[N]) * XNXN - Kprime, Isigma))
# Z = np.dot(Xinput[N], C[N]) * XNXN - Kprime - np.dot(np.transpose(Xinput), np.dot(MSigma, Xinput))
# Z = np.dot(np.transpose(Xinput), np.dot(MSigma - MSigmazeroed, Xinput))
# Z = np.dot(np.dot(Isigma, np.dot(np.transpose(Xinput), np.dot(MSigmazeroed, Xinput))), Isigma)
# Z = Z[OFFSET:N, OFFSET:N]
# Z = M[(OFFSET + 1):(N + 1), (OFFSET + 1):(N + 1)]
# Z = Mchanged[OFFSET:N, OFFSET:N]
# Msigma = Isigma @ Isigma @ M @ Isigma @ Isigma
# Div = np.multiply(Xsigma + Xsigma.T, Xsigma + Xsigma.T)
# Z = (np.log(Msigma / m) / Div)[OFFSET:N, OFFSET:N]
# Z = ((1.0 - ((Msigma / m - 1.0))/Div)/Div)[OFFSET:N, OFFSET:N]
# Z = ((1.0 - np.log(Msigma / m)/Div)/np.log(Sigma))[OFFSET:N, OFFSET:N]
# Z = (np.log(Sigma)/Div)[OFFSET:N, OFFSET:N]
# Z = -np.log(Sigma)[OFFSET:N, OFFSET:N]
# Z = Msigma[OFFSET:N, OFFSET:N]
# Z = np.linalg.inv(M[OFFSET:N, OFFSET:N])
# Z = (Xinput[(1 + OFFSET):, OFFSET:N] @ np.linalg.inv(Kprime[OFFSET:N, OFFSET:N]) @ np.transpose(Xinput[(1 + OFFSET):, OFFSET:N]))
# print(np.linalg.inv(Kprime[:N, :N]))
# Z = -2.0*np.log(Sigma)[OFFSET:N, OFFSET:N]
# Z = (Isigma2 @ M @ Isigma2)[(OFFSET + 1):(N + 1), (OFFSET + 1):(N + 1)]
# Z = np.log(np.abs(MSigmazeroed[(OFFSET + 1):, (OFFSET + 1):]))
# Z = np.log(np.abs(M[(OFFSET + 1):, (OFFSET + 1):]))
# Z = np.log(np.abs(MSigma[(OFFSET + 1):, (OFFSET + 1):]))
# Z = np.log(np.abs(MSigma[OFFSET:N, OFFSET:N]))

# W = np.zeros((N + 1 - OFFSET, N + 1 - OFFSET), dtype='double')
# for i in range(OFFSET, N + 1):
# 	for j in range(OFFSET, N + 1):
# 		if i == N or j == N:
# 			W[i - OFFSET][j - OFFSET] = -M[i][j]*Sigma[i][j]/m
# 		else:


# OFFSET <= J < N
# J = OFFSET
# In = np.zeros((N - J, N - J), dtype='double')
# for i in range(J + 1, N):
# 	In[i - (J + 1)][i - (J + 1)] = np.sqrt(-1.0 / (INPUT[INDEX + 1][i] * M[i][i]))
# In[N - J - 1][N - J - 1] = np.sqrt(1.0 / np.dot(Xinput[i], C[i]))

# A = In @ MSigma[(J + 1):(N + 1), (J + 1):(N + 1)] @ In
# print(A)
# (evals, evecs) = np.linalg.eig(A)
# print(evals)
# print(evecs)

# print(np.diag(Z))
# print(np.dot(Xinput[N], C[N]))

# for i in range(N - OFFSET):
# 	Z[i][i] = 0.0

# print(Z)
# print(Z @ np.ones((N - OFFSET), dtype='double'))

# print(np.linalg.det(Xinput[4:N:(N//4-1), 0:4]))

# Z = (Isigma2 @ M @ Isigma2)[OFFSET:N, OFFSET:N]
# Z = np.linalg.inv(M[OFFSET:N, OFFSET:N])

# q = np.zeros((N - OFFSET), dtype='double')
# w = np.zeros((N - OFFSET), dtype='double')
# for i in range(N - OFFSET):
# 	q[i] = C[N][i + OFFSET] / Xinput[N][i + OFFSET]
# for i in range(N - OFFSET - 1):
# 	# w[i] = 1.0/(q[i + 1] - q[i])
# 	# w[i] = (1.0/(q[i + 1] - q[i]))/np.sqrt(INPUT[INDEX + 1][i + OFFSET] * INPUT[INDEX + 1][i + OFFSET + 1])
# 	w[i] = 1.0/(C[N][OFFSET + i + 1] * Xinput[N][OFFSET + i] - C[N][OFFSET + i] * Xinput[N][OFFSET + i + 1])
# 	# w[i] = (1.0/(C[N][OFFSET + i + 1] * Xinput[N][OFFSET + i] - C[N][OFFSET + i] * Xinput[N][OFFSET + i + 1]))/(INPUT[INDEX + 1][i + OFFSET] * INPUT[INDEX + 1][i + OFFSET + 1])
# print(w)
# print(np.array(INPUT[INDEX + 1])[OFFSET:N])

# for k in range(OFFSET + 1, N - 1):
# 	v = (C[N][k] - C[N][k - 1]) * (Xinput[N][k] - Xinput[N][k + 1]) - (C[N][k] - C[N][k + 1]) * (Xinput[N][k] - Xinput[N][k - 1])
# 	v /= (C[N][k] * Xinput[N][k + 1] - C[N][k + 1] * Xinput[N][k]) * (C[N][k - 1] * Xinput[N][k] - C[N][k] * Xinput[N][k - 1])
# 	print(v)

# Z = Xinput[OFFSET:N, OFFSET:N] @ np.diag(1.0 / Xinput[N, OFFSET:N])
# Z = C[OFFSET:N, OFFSET:N] @ np.diag(1.0 / Xinput[N, OFFSET:N])
# Z = C[OFFSET:N, OFFSET:N] @ np.diag(1.0 / C[N, OFFSET:N])

# for i in range(OFFSET, N):
# 	print(-M[i][i] * INPUT[INDEX + 1][i]**2, C[N][i] * Xinput[N][i] * INPUT[INDEX + 1][i])
# v = 0.0
# f = C[N][N-1]/Xinput[N][N-1]
# print("m:", m)
# for i in range(OFFSET, N):
# 	print(C[N][i]/Xinput[N][i])
# 	# print(v - Csigma[i][i] * Csigma[N][i] * Xsigma[N][i])
# 	# v = Csigma[i][i] * Csigma[N][i] * Xsigma[N][i]
# 	# print(Csigma[N][i] * Csigma[i + 1][i] / Iexp[i][i]**2 * Xsigma[N][i])
# 	# print(np.dot(Xinput[N, (i+1):], C[N, (i+1):])/(f - Csigma[N][i]/Xsigma[N][i]), '\t', f - C[N][i]/Xinput[N][i])
# 	# print(Xsigma[N][(i+1)]*Csigma[N][(i+1)]*(Csigma[N][i + 1]/Xsigma[N][i + 1] - Csigma[N][i]/Xsigma[N][i])/(f - Csigma[N][i]/Xsigma[N][i]), '\t', f - C[N][i]/Xinput[N][i])
# 	# print(Xsigma[N][(i+1)]*Csigma[N][(i+1)]/(Csigma[N][i+1]*Xsigma[N][i] - Csigma[N][i]*Xsigma[N][i+1]))
# print(np.dot(C[N], Xinput[N]) + 2.0 * np.exp(-np.dot(Xinput[N], Xinput[N]) / 2.0))


Cinv = np.linalg.inv(C)
ID = np.eye(N + 1, N + 1, dtype='double')
for i in range(OFFSET, N):
	for k in range(OFFSET, i + 1):
		Z[i - OFFSET][k - OFFSET] = C[N, (i + 1):] @ Cinv[(i + 1):, (k + 1):] @ (-INPUT[INDEX + 1][k + 1] * ID[(k + 1):, k + 1] + C[(k + 1):, (k + 1):] @ (C[k + 1, (k + 1):]).T)
		Z[i - OFFSET][k - OFFSET] = Z[i - OFFSET][k - OFFSET]**2 * (Xinput[N][i] / C[N][i] - Xinput[N][i + 1] / C[N][i + 1])
		# if i == N - 1:
		# 	Z[i - OFFSET][k - OFFSET] = 0.0
		if k < i:
			Z[i - OFFSET][k - OFFSET] += Z[i - OFFSET - 1][k - OFFSET]
		# if i == N - 1:
		# 	Z[i - OFFSET][k - OFFSET] = 0.0

for j in range(OFFSET, N):
	# print(Xinput[N][j] / C[N][j])
	print(Xinput[N][j] / C[N][j] - Xinput[N][j + 1] / C[N][j + 1])

print("TOTAL:")
print(Xinput[N][OFFSET] / C[N][OFFSET])

# ==============================================================================================================================

# Z = np.array([[Cinv[i + OFFSET][j + OFFSET] * np.sqrt(INPUT[INDEX + 1][i + OFFSET]) for j in range(MAX - OFFSET)] for i in range(MAX - OFFSET)])
# Z = np.array([[C[i + 1 + OFFSET][j + OFFSET] * np.sqrt(INPUT[INDEX + 1][j + OFFSET])  for j in range(MAX - OFFSET)] for i in range(MAX - OFFSET)])
# Z = np.array([[INPUT[INDEX][i + 1 + OFFSET][j + OFFSET] * np.sqrt(INPUT[INDEX + 1][j + OFFSET]) for j in range(MAX - OFFSET)] for i in range(MAX - OFFSET)])
# Z = np.array([[np.log(INPUT[INDEX][i + 1 + OFFSET][j + OFFSET] * np.sqrt(INPUT[INDEX + 1][j + OFFSET])) for j in range(MAX - OFFSET)] for i in range(MAX - OFFSET)])

# for i in range(N - OFFSET):
# 	for j in range(N - OFFSET):
# 		if j >= i - 11 or i >= N - 11:
# 			Z[i][j] = 0.0

# for i in range(N - OFFSET):
# 	if i == 0:
# 		Z[i][i] = Z[i + 1][i]
# 	else:
# 		Z[i][i] = Z[i][i - 1]

# for i in range(N - OFFSET):
# 	for j in range(N - OFFSET):
# 		if abs(i - j) <= 1 or i == 0 or j == 0:
# 			Z[i][j] = 0.0
# 		# if i == 0 or j == 0:
# 		# 	Z[i][j] = 0.0
# 		if abs(i - j) <= 1:
# 			Z[i][j] = 0.0

# for i in range(N - OFFSET):
# 	for j in range(N - OFFSET):
# 		if abs(i - j) <= 0:
# 			if i > 0 and j > 0 and i < N - OFFSET - 1 and j < N - OFFSET - 1:
# 				Z[i][j] = (Z[i - 1][j] + Z[i][j - 1] + Z[i + 1][j] + Z[i][j + 1]) / 4.0
# 			else:
# 				Z[i][j] = 0.0

# for i in range(MAX - OFFSET):
# 	if i == 0:
# 		Z[i][i] = Z[i + 1][0]
# 	else:
# 		Z[i][i] = Z[i][i - 1]

# for j in range(MAX - OFFSET):
# 	c = Z[j][j]
# 	for i in range(j, MAX - OFFSET):
# 		Z[i][j] /= c
# 		# Z[i][j] = np.log(Z[i][j]/c)

# for i in range(1, MAX - OFFSET):
# 	c = Y[i][i]
# 	for j in range(i + 1):
# 		Y[i][j] = (Y[i][j] * 1.0) / c

for i in range(0, MAX - OFFSET):
	for j in range(i, MAX - OFFSET):
		X[i][j] = X[j][j]
		Y[i][j] = Y[j][j]
		Z[i][j] = Z[j][j]

ax.set_xlabel('used up precision')
ax.set_ylabel('coordinate index')
ax.set_zlabel('X')
# Z = np.array(Z)
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)


plt.show()
