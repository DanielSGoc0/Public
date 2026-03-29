# This program is useful if we don't search for optimal nodes with small error (like 10^(-4)).
# Instead, it is useful if the error is relatively big, say >= 1.

import numpy as np
import scipy
import math

S2 = math.sqrt(2.0)
LS2 = math.log(2.0)/2.0

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
	global S2, LS2
	P = S2*T + math.log((1 + math.exp(-S2*T)) / 2.0)
	STEPS = math.ceil(P/h)

	X = np.array(X, dtype='double')
	N = X.shape[0]

	M = np.zeros((N, N), dtype='double')
	for i in range(N):
		for j in range(N):
			M[i][j] = math.exp(np.dot(X[i], X[j]))
	V = np.linalg.cholesky(M)

	ANS = 2.0*P - 2.0*math.log(2.0 - math.exp(-P))

	# find sequence of functions
	gs = [lambda p: S2 * (1 - math.exp(-p))]
	
	for k in range(1, N):
		a = scipy.linalg.solve_triangular(np.transpose(V), X[k], lower=False)
		gk = lambda p: math.exp(LS2 - p + sum([a[i] * gs[i](p) for i in range(k)]))
		gs.append(create_callable_integral(gk, 0.0, P, STEPS))

	g = make_vector(gs, P, N)
	x = scipy.linalg.solve_triangular(V, g, lower=True)
	ANS -= np.dot(x, x)

	print(ANS)

	return ANS


# =========================================================================================

def initial_guess(N, T, h):
	global S2, LS2
	STEPS = math.ceil(T/h)

	X = np.zeros((N, N), dtype='double')
	gs = [lambda t: S2 * math.tanh(S2 * t / 2.0)]

	for k in range(1, N):

		M = np.zeros((k, k), dtype='double')
		for i in range(k):
			for j in range(k):
				M[i][j] = math.exp(np.dot(X[i], X[j]))
		V = np.linalg.cholesky(M)

		g = make_vector(gs, 1.0 * T * k/N, k)
		b = scipy.linalg.solve_triangular(V, g, lower=True)
		for i in range(k):
			X[k][i] = b[i]

		if k == N - 1:
			break

		a = scipy.linalg.solve_triangular(np.transpose(V), b, lower=False)
		gk = lambda t: math.cosh(S2 * t / 2.0)**(-2) * math.exp(sum([a[i] * gs[i](t) for i in range(k)]))
		gs.append(create_callable_integral(gk, 0.0, T, STEPS))

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

minimize(3, 1, 0.01)
