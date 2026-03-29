# A practical (and naive) implementation of Minimizer resulting from TrueGradient

import numpy as np

D = 4
A = np.zeros((D, D), dtype='double')


def coeffs(n):
	# return 1/np.log(n + 2)
	# return 1.0
	# return 8.0/(n + 6)
	return np.sqrt(8.0/(n + 8))


def f(x):
	global A
	# v = np.dot(A, x)
	# return np.dot(v, v)
	# return np.sin(x[0])
	return -((1 - x[0])**2 + 100*(x[1] - x[0]**2)**2)
	# return -x[0]**2

def df(x):
	global A
	# return 2*np.dot(np.transpose(A), np.dot(A, x))
	# return np.array([np.cos(x[0])])
	return -np.array([-2*(1 - x[0]) - 400*x[0]*(x[1] - x[0]**2), 200*(x[1] - x[0]**2)])
	# return np.array([-2*x[0]])


# minimize the function f
# df is the derivative of f.
# x0 is the starting point
# N is the amoint of step (must be >= 1)
def maximize(f, df, x0, h, N):
	x0 = np.array(x0)
	d = len(x0)

	a = [coeffs(i) for i in range(N + 1)]
	w = [np.exp(-a[i]**2 / 2.0) for i in range(N + 1)]
	S = [1.0 for i in range(N + 1)]
	for i in range(1, N + 1):
		S[i] = np.sqrt(1 - w[i - 1]**2)

	x = np.zeros((N + 1, d), dtype='double')
	F = [0.0 for i in range(N + 1)]
	DF = np.zeros((N + 1, d), dtype='double')
	e = np.zeros((N + 1, d), dtype='double')

	x[0] = x0
	F[0] = 0.0
	DF[0] = df(x0)
	e[0] = DF[0]
	print(x[0], f(x[0]))

	x[1] = x[0] + h * a[0] * DF[0]
	F[1] = w[0] * a[0]
	DF[1] = df(x[1])
	e[1] = 1.0/S[1] * (DF[1] + a[0] * F[1] * e[0] - w[0] * DF[0])
	print(x[1], f(x[1]))

	for n in range(1, N):
		x[n + 1] = x[n] + h * a[n] * e[n]
		# F[n + 1] = w[n] * (a[n] * S[n] + F[n])
		F[n + 1] = f(x[n + 1])
		DF[n + 1] = df(x[n + 1])
		e[n + 1] = 1/S[n + 1] * (DF[n + 1] - 1.0/S[n] * w[n - 1]**2 * w[n] * a[n - 1] * a[n] * e[n - 1] + a[n] * F[n + 1] * e[n] - w[n] * DF[n])

		# print(x[n + 1], F[n + 1], DF[n + 1], e[n + 1], f(x[n + 1]))
		print(x[n + 1], f(x[n + 1]))



N = 100

# maximize(f, df, [0.0], 1.0, N)
# maximize(f, df, [1.0], 0.1, N)
maximize(f, df, [0.0, 0.0], 0.01, N)
# create_random_matrix()
# print(A)
# A = np.identity(N)
# for i in range(N):
# 	A[i][i] = i + 1
# current_x = np.ones(N)
# maximize()
