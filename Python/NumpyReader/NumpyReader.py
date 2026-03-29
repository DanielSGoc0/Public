# A file used to read data from file numpy_input.txt

import numpy as np
from pathlib import Path
import scipy

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:5.8f}'.format}, linewidth=200)

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

X = np.array(read_from_file("numpy_input.txt")[0])
# print(X)

# ==================================================================================

# def kronecker(i, j):
# 	if i == j:
# 		return 1.0
# 	else:
# 		return 0.0
	
# def changer(x, y):
# 	# return np.exp((x-1)**2)
# 	# return (x-1.3)**2
# 	# return 1.0
# 	# return (4.0 - x**2)/4.0
# 	# return (x-2)**2/4.0
# 	return 1.0
# 	# return y
# 	# return (4.0 - x**2)/8.0 + y
	
# N = X.shape[0] - 1
N = X.shape[0]
err2 = 1.0

# Geo = np.zeros((N + 1, N + 1), dtype='double')
# for i in range(N + 1):
# 	for j in range(N + 1):
# 		Geo[i][j] = np.dot(X[i] - X[j], X[i] - X[j])
# print(Geo)


# Sigma = np.zeros((N + 1, N + 1), dtype='double')
# for i in range(N + 1):
# 	for j in range(N + 1):
# 		Sigma[i][j] = np.exp(-np.dot(X[i] - X[j], X[i] - X[j])/2.0)
# 	Sigma[i][i] += err2
# print(scipy.linalg.cholesky(Sigma, lower=True))

# X = np.array(X)
# X2 = np.zeros((N + 1), dtype='double')
# F = np.zeros((N + 1), dtype='double')
# DF = np.zeros((N + 1, N + 1), dtype='double')
# DF2 = np.zeros((N + 1), dtype='double')
# Sigma = np.zeros((N + 1, N + 1), dtype='double')
# Inverses_before = [np.zeros((N + 1, N + 1), dtype='double') for i in range(N + 2)]
# Inverses_after = [np.zeros((N + 1, N + 1), dtype='double') for i in range(N + 2)]

# # initial values
# DF2[0] = 1.0
# DF[0][0] = np.sqrt(1.0 + changer(F[0], DF2[0]) * err2)
# Sigma[0][0] = 1.0


# for n in range(1, N + 1):
# 	X2[n] = np.dot(X[n], X[n])

# 	for i in range(n + 1):
# 		Sigma[n][i] = np.exp(-(X2[n] + X2[i]) / 2.0 + np.dot(X[i], X[n]))
# 		Sigma[i][n] = Sigma[n][i]

# 	# Next, calculate the magnitude^2 of orthogonal vector.
# 	v = np.array(Sigma[n, :n])
# 	Sigma_partial = np.array(Sigma[:n, :n])

# 	for i in range(n):
# 		Sigma_partial[i][i] += err2*changer(F[i], DF2[i])
# 	c = 1 - np.dot(v, scipy.linalg.solve(Sigma_partial, v, assume_a='pos'))
# 	# print("vector:")
# 	print(scipy.linalg.solve(Sigma_partial, v, assume_a='pos'))

# 	INV = np.linalg.inv(Sigma_partial)
# 	for i in range(N + 1):
# 		for j in range(N + 1):
# 			if i < n and j < n:
# 				Inverses_after[n][i][j] = INV[i][j]


# 	# And then calculate the projection of DF[n] onto first n vectors.
# 	M = np.zeros((n*n, n*n), dtype='double')
# 	for a in range(n):
# 		for b in range(n):
# 			for i in range(n):
# 				for j in range(n):
# 					M[a*n + i][b*n + j] = Sigma[a][b] * (kronecker(i, j) - (X[a][i] - X[b][i]) * (X[a][j] - X[b][j]))
# 	for a in range(n):
# 		for i in range(n):
# 			M[a*n + i][a*n + i] += err2*changer(F[a], DF2[a])

# 	d = np.zeros(n*n, dtype='double')
# 	for a in range(n):
# 		for i in range(n):
# 			d[a*n + i] = DF[a][i]
# 	r = scipy.linalg.solve(M, d, assume_a='pos')
# 	print("r = ")
# 	print(r)
# 	# print("M inverse:")
# 	# print(np.linalg.inv(M))

# 	L = np.zeros((n, n*n), dtype='double')
# 	for b in range(n):
# 		for i in range(n):
# 			for j in range(n):
# 				L[i][b*n + j] = Sigma[n][b] * (kronecker(i, j) - (X[n][i] - X[b][i]) * (X[n][j] - X[b][j]))

# 	mean = np.dot(L, r)


# 	# As well as the mean of F(Xn).
# 	l = np.zeros(n*n, dtype='double')
# 	for b in range(n):
# 		for j in range(n):
# 			l[b*n + j] = Sigma[n][b] * (X[n][j] - X[b][j])

# 	F[n] = np.dot(l, r)


# 	# Finally, update the vector DF[n] accordingly.
# 	for i in range(n):
# 		DF[n][i] = mean[i]
# 	DF2[n] = c + np.dot(DF[n], DF[n])
# 	DF[n][n] = np.sqrt(c + err2*changer(F[n], DF2[n]))

# # Inverses_before[N + 1] = np.linalg.inv(Sigma)
# # print("INV before:")
# # print(Inverses_before[N + 1])

# for i in range(N + 1):
# 	Sigma[i][i] += err2*changer(F[i], DF2[i])
# Inverses_after[N + 1] = np.linalg.inv(Sigma)
# # print("INV after:")
# # print(Inverses_after[N + 1])


# print("==================   RESULTS:   ==================")
# print("F = ")
# print(F)
# print("X = ")
# print(X)
# print("DF = ")
# print(DF)
# # print("Sigma inverse = ")
# # print(np.linalg.inv(Sigma))
# # for i in range(N + 1):
# # 	print(np.dot(DF[i], DF[i]))
# # print("other...")
# # for i in range(1, N):
# # 	print(DF[i][i]**2, 1 - np.exp(-A[i - 1]**2))
# # for i in range(2, N + 1):
# # 	print(DF[i - 1][i - 2], - A[i - 1] * A[i - 2] * np.exp(-A[i - 2]**2)/np.sqrt(1 - np.exp(-A[i - 2]**2)))
# # for i in range(2, N):
# # 	print(A[i - 1] * (A[i - 1] * np.sqrt(1 - np.exp(-A[i - 2]**2)) - A[i]  * np.exp(-A[i - 1]**2)/np.sqrt(1 - np.exp(-A[i - 1]**2))))
# # print("Formula for T:")
# # for k in range(1, N):
# # 	SUM = 0.0
# # 	for i in range(k):
# # 		SUM += Sigma[k][i] * DF[i][i] * A[i]
# # 	print(DF[k][k - 1], np.exp(-A[k - 1]**2 / 2) * DF[k - 1][k - 1] - A[k - 1] * SUM)
# # print("General formula:")
# # for k in range(2, N):
# # 	print(DF[k - 1][k - 1]/A[k - 1] + A[k] * np.exp(-A[k - 1]**2 / 2) / DF[k][k], A[k - 1] / DF[k - 1][k - 1] + np.exp(-A[k - 2]**2 / 2)/A[k - 2] * DF[k - 2][k - 2])
# # print("other....")
# # for n in range(N):
# # 	print(F[n + 1], np.exp(-A[n]**2 / 2) * (DF[n][n] * A[n] + F[n]))

# # print("DFs subtracted:")
# # ei = np.zeros((N + 1))
# # ei[1] = 1.0
# # ei_1 = np.zeros((N + 1))
# # ei_1[0] = 1.0
# # print(DF[1] - np.exp(-A[0]**2 / 2.0) * DF[0] - DF[1][1] * ei + A[0] * F[1] * ei_1)
# # for i in range(2, N + 1):
# # 	ei = np.zeros((N + 1))
# # 	ei[i] = 1.0
# # 	ei_1 = np.zeros((N + 1))
# # 	ei_1[i-1] = 1.0
# # 	ei_2 = np.zeros((N + 1))
# # 	ei_2[i-2] = 1.0
# # 	print(DF[i] - np.exp(-A[i - 1]**2 / 2.0) * DF[i - 1] - DF[i][i] * ei + A[i - 1] * F[i] * ei_1 - 1/DF[i - 1][i - 1]*np.exp(-A[i-2]**2 - A[i-1]**2 / 2.0)*A[i - 1]*A[i - 2]*ei_2)

# # for i in range(1, N + 1):
# # 	print(np.sqrt(-np.log(DF[i][i] / DF[i-1][i-1]) * 2))

# print()
# print(np.dot(DF, np.transpose(DF)))

# for i in range(N + 1):
# 	DF[i][i] = 0.0

# print(np.dot(DF, np.transpose(DF)))



w = np.zeros((N + 1), dtype='double')
F = np.zeros((N + 1), dtype='double')

# initial values
w[0] = 1.0
F[0] = 0.0

for n in range(N):
	# calculate w[n + 1] and F[n + 1]
	w[n + 1] = 1 - np.exp(-X[n]**2) * (1 - err2 * w[n] / (err2 + w[n]))
	F[n + 1] = np.exp(-X[n]**2 / 2.0) * (F[n] + X[n] * w[n] / np.sqrt(err2 + w[n]))

for n in range(0, N + 1):
	# print(-np.log(w[n]) / np.log(n), -np.log(2 - F[n]) / np.log(n))
	print(X[n], w[n], F[n])
