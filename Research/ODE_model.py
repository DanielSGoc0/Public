# This program computes the error in a certain mathematical model.
# The details are bery technical.
# The project can help create a new ML optimizer.

import numpy as np
import scipy
import math

# This variable is the measurement error.
err2 = 0.00

def get_v(N, X):
	v = []
	for i in range(N):
		for j in range(i):
			v.append(X[i][j])
	return v


def get_X(N, v):
	X = np.zeros((N, N), dtype='double')
	k = 0
	for i in range(N):
		for j in range(i):
			X[i][j] = v[k]
			k += 1
	return X


def make_vector(funs, x, k):
	res = np.zeros((k), dtype='double')
	for i in range(k):
		res[i] = funs[i](x)
	return res


# ===============================================================================

def create_callable_integral(f, a, b, nodes):
	xs = []
	ys = []
	dys = []

	xs.append(0)
	ys.append(0)
	dys.append(f(0))
	for i in range(1, nodes + 1):
		x = a + (b - a) * i / nodes
		xs.append(x)
		ys.append(ys[-1] + scipy.integrate.fixed_quad(lambda v: [f(x) for x in v], xs[-2], xs[-1], n=5)[0])
		dys.append(f(x))

	return scipy.interpolate.CubicHermiteSpline(xs, ys, dys)


def calc(X, T, h):
	global err2
	P = 2*T + math.log((1 + math.exp(-2*T))/2.0)
	STEPS = math.ceil(P/h)

	X = np.array(X, dtype='double')
	N = X.shape[0]

	M = np.zeros((N, N), dtype='double')
	for i in range(N):
		for j in range(N):
			M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]))
		M[i][i] += err2
	print(X)
	V = np.linalg.cholesky(M)

	ANS = P - math.log(2.0 - math.exp(-P))

	# find sequence of functions
	gs = [lambda p: 1 - math.exp(-p)]
	L2 = [np.dot(X[i], X[i]) for i in range(N)]

	for k in range(1, N):
		a = scipy.linalg.solve_triangular(np.transpose(V), X[k], lower=False)
		gk = lambda p: math.exp(-p + 2*sum([a[i] * gs[i](p) for i in range(k)]) - L2[k])
		gs.append(create_callable_integral(gk, 0.0, P, STEPS))

	g = make_vector(gs, P, N)
	x = scipy.linalg.solve_triangular(V, g, lower=True)
	ANS -= np.dot(x, x)

	print(ANS)

	return ANS


# =========================================================================================

def initial_guess(N, T, h):
	global err2
	P = 2*T + math.log((1 + math.exp(-2*T))/2.0)
	STEPS = math.ceil(P/h)

	X = np.zeros((N, N), dtype='double')
	gs = [lambda p: 1 - math.exp(-p)]

	for k in range(1, N):

		M = np.zeros((k, k), dtype='double')
		for i in range(k):
			for j in range(k):
				M[i][j] = math.exp(-np.dot(X[i] - X[j], X[i] - X[j]))
			M[i][i] += err2
		V = np.linalg.cholesky(M)

		g = make_vector(gs, 1.0 * P * k/N, k)
		b = scipy.linalg.solve_triangular(V, g, lower=True)
		for i in range(k):
			X[k][i] = b[i]

		if k == N - 1:
			break

		L2 = np.dot(b, b)
		a = scipy.linalg.solve_triangular(np.transpose(V), b, lower=False)
		gk = lambda p: math.exp(-p + 2*sum([a[i] * gs[i](p) for i in range(k)]) - L2)
		gs.append(create_callable_integral(gk, 0.0, P, STEPS))

	return X



def minimize(N, T, h):
	v0 = initial_guess(N, T, h)
	print(v0)

	fun = lambda v: calc(get_X(N, v), T, h)

	RES = scipy.optimize.minimize(fun, get_v(N, v0))
	X = get_X(N, RES.x)
	print(RES)
	for i in range(N - 1):
		print(np.dot(X[i + 1] - X[i], X[i + 1] - X[i]))


# print(calc([[0.0, 0.0, 0.0], [0.4385, 0.0, 0.0], [0.6541, 0.2786, 0.0]], 10, 0.01))

minimize(10, 10, 0.01)
