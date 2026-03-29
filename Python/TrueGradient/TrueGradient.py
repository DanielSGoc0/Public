# The covariance function is e^(-||x - y||^2 / 2)
# No modifications. This problem has been essentially solved.
import scipy
import numpy as np
import math
import random


# =========================================================================================

def kronecker(i, j):
	if i == j:
		return 1.0
	else:
		return 0.0

def get_v(N, A):
	v = []
	for i in range(N):
		for j in range(i + 1):
			v.append(A[i][j])
	return v

def get_A(N, v):
	A = np.zeros((N, N), dtype='double')
	k = 0
	for i in range(N):
		for j in range(i + 1):
			A[i][j] = v[k]
			k += 1
	return A


# =========================================================================================

# A is lower triangular (not strictly)
def calc(A):
	print("A = ", A)
	A = np.array(A, dtype='double')
	N = A.shape[0]

	# X describes the geometry of points Xn and is strictly lower triangular
	# X2 is a list of magnitudes squared of vectors in X.
	# F describes the means of F(Xn)/sqrt(dim).
	# DF describes the geometry of gradients dF(Xn) and is lower triangular
	# Sigma_ij = e^(-||Xi - Xj||^2 / 2)
	X = np.zeros((N + 1, N + 1), dtype='double')
	X2 = np.zeros((N + 1), dtype='double')
	F = np.zeros((N + 1), dtype='double')
	DF = np.zeros((N + 1, N + 1), dtype='double')
	Sigma = np.zeros((N + 1, N + 1), dtype='double')

	# initial values
	DF[0][0] = 1.0
	Sigma[0][0] = 1.0


	for n in range(1, N + 1):
		# At n-th step we know (X[0], ... , X[n - 1]), (DF[0], ... , DF[n - 1])
		# We wish to learn the points X[n] and DF[n].

		# Firstly, we update the geometry of X[n] and Sigma:
		# X[n] = A[n - 1][0] * DF[0] + ... + A[n - 1][n - 1] * DF[n - 1]
		for i in range(n):
			X[n] += A[n - 1][i] * DF[i]

		X2[n] = np.dot(X[n], X[n])

		for i in range(n + 1):
			Sigma[n][i] = math.exp(-(X2[n] + X2[i]) / 2.0 + np.dot(X[i], X[n]))
			Sigma[i][n] = Sigma[n][i]


		# Next, calculate the magnitude^2 of orthogonal vector.
		v = np.array(Sigma[n, :n])
		Sigma_partial = np.array(Sigma[:n, :n])
		c = 1 - np.dot(v, scipy.linalg.solve(Sigma_partial, v, assume_a='pos'))


		# And then calculate the projection of DF[n] onto first n vectors.
		M = np.zeros((n*n, n*n), dtype='double')
		for a in range(n):
			for b in range(n):
				for i in range(n):
					for j in range(n):
						M[a*n + i][b*n + j] = Sigma[a][b] * (kronecker(i, j) - (X[a][i] - X[b][i]) * (X[a][j] - X[b][j]))

		d = np.zeros(n*n, dtype='double')
		for a in range(n):
			for i in range(n):
				d[a*n + i] = DF[a][i]
		r = scipy.linalg.solve(M, d, assume_a='pos')

		L = np.zeros((n, n*n), dtype='double')
		for b in range(n):
			for i in range(n):
				for j in range(n):
					L[i][b*n + j] = Sigma[n][b] * (kronecker(i, j) - (X[n][i] - X[b][i]) * (X[n][j] - X[b][j]))

		mean = np.dot(L, r)


		# As well as the mean of F(Xn).
		l = np.zeros(n*n, dtype='double')
		for b in range(n):
			for j in range(n):
				l[b*n + j] = Sigma[n][b] * (X[n][j] - X[b][j])

		F[n] = np.dot(l, r)


		# Finally, update the vector DF[n] accordingly.
		for i in range(n):
			DF[n][i] = mean[i]
		DF[n][n] = math.sqrt(c)

	print("F = ", F)
	# print("Geo = ", np.dot(DF, np.transpose(DF)))
	# G = [[np.dot(X[i] - X[j], X[i] - X[j]) for j in range(N + 1)] for i in range(N + 1)]
	# print(np.array(G))
	print(X)
	return F[N]


# =========================================================================================

def initial_guess(N):
	A = np.zeros((N, N), dtype='double')

	for i in range(N):
		for j in range(N):
			if j <= i:
				A[i][j] = random.random()
			else: 
				A[i][j] = 0.0

	return A


def minimize(N):
	A0 = initial_guess(N)
	print(A0)

	fun = lambda v: -calc(get_A(N, v))

	RES = scipy.optimize.minimize(fun, get_v(N, A0))
	X = get_A(N, RES.x)
	print(RES)


# minimize(5)
calc([[0, 0, 0, 0], [0.5, 0, 0, 0], [0.3, 0.8, 0, 0], [0.7, -0.2, 0.4, 0]])


# N = 1:
# F =  	[0.        	0.60653066]
# A =  [[1.00000002]]


# N = 2:
# F =  	[0.         0.58065305	0.93619482	]
# A =  [[1.2155948  0.        ]
# 	   [1.40329494  0.8226424 ]]


# N = 3:
# F =  	[0.        	0.5680322	0.90274564	1.14070388	]
# A =  [[1.26612912	0.        	0.        ]
#		[1.55340746	1.06174748	0.        ]
#		[1.61800998	1.30051765	0.74387267]]


# N = 4:
# F = 	[0.        	0.56422545	0.88126505	1.10755178	1.27966353	]
# A =  [[1.27992269 0.        	0.        	0.        ]
#		[1.59816859	1.1311675 	0.        	0.        ]
#		[1.70950855	1.52693572	0.99038043	0.        ]
#		[1.7383227 	1.6293565 	1.24667005	0.6974091 ]]


# N = 5:
# F =  	[0.         0.563027   	0.87354612 	1.08169155 	1.24922344 	1.3801858	]
# A =  [[1.28415812	0.        	0.        	0.        	0.       	]
# 		[1.61226956	1.1529521 	0.        	0.        	0.       	]
# 		[1.74125361	1.60620431	1.07256782	0.        	0.       	]
#		[1.79528739 1.79607411	1.52186129	0.94511039	0.       	]
#		[1.81021919	1.84853344	1.64599064 	1.20619861 	0.66629011]	]


# N = 6:
# F =  	[0.         0.56260407 	0.8708163  	1.07104311 	1.22173565 	1.35295979 	1.45626668]
# A =  [[1.28564136 0.         	0.         	0.         	0.         	0.        ]
# 		[1.61715214 1.16035065 	0.         	0.         	0.         	0.        ]
# 		[1.7524979  1.63404082 	1.10111951 	0.         	0.         	0.        ]
# 		[1.81750803 1.8615544  	1.62999237 	1.03609886 	0.         	0.        ]
# 		[1.8474298  1.9662627  	1.87339717 	1.51295367 	0.91263355 	0.        ]
# 		[1.85597431 1.99616156 	1.9428991  	1.64911759 	1.1732276  	0.64381518]]


# N = 7:
# F =  	[0.         0.56246364 	0.86978365 	1.06690914 	1.20906793 	1.32535419 	1.43212533 	1.51584616]
# A =  [[1.28613258 0.         	0.         	0.         	0.         	0. 			0.        ]
#		[1.61888158 1.16313576 	0.         	0.         	0.         	0.			0.        ]
#		[1.75657298 1.64448363 	1.11160906 	0.         	0.         	0.			0.        ]
#		[1.82592839 1.88694125 	1.6714144  	1.07037399 	0.         	0.			0.        ]
# 		[1.86315155 2.01706198 	1.97179899 	1.64473491 	1.00984557 	0.			0.        ]
#		[1.8812394  2.08030168 	2.11777976 	1.92387509 	1.50066459 	0.88763023	0.        ]
#		[1.8864939  2.09867525 	2.16019039 	2.00497563 	1.64325905 	1.14547449	0.62679059]]


# N = 8:
# F =  	[0.         0.5622378  	0.86937022 	1.06543403 	1.20379354 	1.3113478	1.40528003 	1.49450302 	1.56376149]
# A =  [[1.28692121 0.         	0.         	0.         	0.         	0.			0.         	0.        ]
#		[1.62065426 1.16403535 	0.         	0.         	0.         	0.			0.         	0.        ]
# 		[1.75985964 1.64790201 	1.11520859 	0.         	0.         	0.			0.         	0.        ]
#		[1.83153612 1.89604844 	1.68693122 	1.08399993 	0.         	0.			0.         	0.        ]
#		[1.87235233 2.03678749 	2.01115922 	1.6993892  	1.04911731 	0.			0.         	0.        ]
# 		[1.89585081 2.11754874 	2.19721476 	2.05277198 	1.65178856 	0.98896518	0.         	0.        ]
#		[1.90767823 2.15812285 	2.29068789 	2.23037658 	1.95471317 	1.48584237	0.86749836 	0.        ]
# 		[1.91115452 2.17003229 	2.31812182 	2.28250486 	2.04361524 	1.63160524	1.12182794 	0.61358788]]
